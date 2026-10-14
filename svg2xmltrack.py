#!/usr/bin/env python3

import argparse
import os
from svgpathtools import svg2paths
from svgpathtools.path import Line, Path
import numpy as np
from scipy.optimize import least_squares
from lxml import etree
import xmlformatter


def _curve_nb_steps(R_s, R_e, arc, step_length):
    """Number of steps TrackGen divides a curve into (radii and step_length in meters)."""
    length = arc * (R_s + R_e) / 2
    return int(length / step_length) + 1


def _simulate_curve_endpoint(R_s, R_e, arc, n_steps):
    """Simulate how TrackGen builds a curve with n_steps steps.
    Returns (x_fwd, y_ctr): chord components along and perpendicular to the start
    heading, where y_ctr > 0 is the side of the curve center.
    """
    length = arc * (R_s + R_e) / 2
    step_length = length / n_steps

    if R_e != R_s and n_steps != 1:
        radii = R_s + (R_e - R_s) / (n_steps - 1) * np.arange(n_steps)
        # TrackGen rescales step_length so the total rotation equals arc exactly
        step_length *= arc / np.sum(step_length / radii)
    else:
        # With a single step TrackGen uses only the start radius and does not rescale
        radii = np.full(n_steps, R_s)

    # Each step is an exact circular arc of radius radii[k] turning step_arcs[k]
    step_arcs = step_length / radii
    end_angles = np.cumsum(step_arcs)
    start_angles = end_angles - step_arcs
    x_fwd = np.sum(radii * (np.sin(end_angles) - np.sin(start_angles)))
    y_ctr = np.sum(radii * (np.cos(start_angles) - np.cos(end_angles)))
    return x_fwd, y_ctr


def calculate_radius_end_radius(angle_radians, segment, prev_segment, scale=1.0, step_length=4.0):
    """Find (radius, end radius) so that TrackGen's curve ends where the SVG segment ends.

    scale is the number of meters per SVG unit and step_length is the profil steps length
    in meters. TrackGen's step count depends on the curve length in meters, so the
    simulation runs in meters; the result is in SVG units.
    Returns None when no single curve can end there (see calculate_biarc).
    """
    arc = abs(angle_radians)
    if arc < 1e-4:
        return None

    # Chord (in meters) in the local frame aligned with the start heading
    chord = (segment.end - segment.start) * _heading(prev_segment).conjugate() * scale
    target = chord.real, chord.imag * np.sign(angle_radians)

    # Start from a constant radius. The step count is an integer that depends on the radii,
    # so solve with it fixed and repeat until stable. The bounds keep the radii positive
    # and finite, since unreachable endpoints push them to infinity.
    radii = np.full(2, abs(chord) / (2 * np.sin(arc / 2)))
    bounds = (1e-6 * abs(chord), 20 * abs(chord) / arc)
    for _ in range(10):
        n_steps = _curve_nb_steps(*radii, arc, step_length)
        radii = least_squares(lambda r: np.subtract(_simulate_curve_endpoint(*r, arc, n_steps), target),
                              np.clip(radii, *bounds), bounds=bounds, x_scale="jac", xtol=1e-15, ftol=1e-15, gtol=1e-15).x
        if _curve_nb_steps(*radii, arc, step_length) == n_steps:
            break

    R_s, R_e = radii
    n_steps = _curve_nb_steps(R_s, R_e, arc, step_length)
    # Unreachable endpoint (e.g. S-shaped SVG curve)
    if np.hypot(*np.subtract(_simulate_curve_endpoint(R_s, R_e, arc, n_steps), target)) > 1e-4 * abs(chord):
        return None
    # A single step turns length / R_s instead of arc, so it only works with a constant radius
    if n_steps == 1 and abs(arc * (R_e - R_s) / (2 * R_s)) > 1e-7:
        return None
    return R_s / scale, R_e / scale

def calculate_biarc(segment, prev_segment, next_segment):
    """Join two straights with two consecutive constant radius curves (a biarc).

    Used when a single TrackGen curve cannot end where the SVG segment ends. Returns two
    (angle_radians, radius, end radius, length) tuples in SVG units, where the sign of the
    angle follows angle_between_segments(). A piece that does not turn is a straight and
    its radius is None. TrackGen builds constant radius curves exactly, whatever its
    step count is.
    """
    def dot(a, b):
        return (a * b.conjugate()).real

    heading_start, heading_end = _heading(prev_segment), _heading(next_segment)
    chord = segment.end - segment.start
    headings = heading_start + heading_end

    # Both curves use tangents of the same length
    c = 1 - dot(heading_start, heading_end)
    if c < 1e-12:
        if dot(chord, heading_end) <= 0:
            raise Exception("Two parallel straights cannot be joined when the second one starts behind the first one")
        tangent = dot(chord, chord) / (4 * dot(chord, heading_end))
    else:
        tangent = (-dot(chord, headings) + np.sqrt(dot(chord, headings) ** 2 + 2 * c * dot(chord, chord))) / (2 * c)

    # Heading at the point where both curves meet
    heading_joint = (segment.end - tangent * heading_end) - (segment.start + tangent * heading_start)
    heading_joint = heading_joint / abs(heading_joint)

    pieces = []
    for heading_in, heading_out in ((heading_start, heading_joint), (heading_joint, heading_end)):
        angle = np.angle(heading_out / heading_in)
        if abs(angle) < 1e-6:
            pieces.append((0, None, None, 2 * tangent))
        else:
            radius = tangent / np.tan(abs(angle) / 2)
            pieces.append((angle, radius, radius, radius * abs(angle)))
    return pieces

def _heading(segment):
    """Unit vector (as a complex number) from the start to the end of a segment."""
    return (segment.end - segment.start) / abs(segment.end - segment.start)

def angle_between_segments(segment1, segment2):
    """Signed angle in radians from segment1 to segment2; positive is a right turn in SVG coordinates."""
    return np.angle(_heading(segment2) / _heading(segment1))

def length_in_meters(value, parameters):
    if parameters["track_length"] == 0:
        return value

    return format_number((value / parameters["svg_length"]) * parameters["track_length"])

def create_xml_file(segments, parameters):
    doctype = """<!DOCTYPE params SYSTEM "../../../src/libs/tgf/params.dtd" [
<!--  general definitions for tracks  -->
<!ENTITY default-surfaces SYSTEM "../../../data/tracks/surfaces.xml">
]>"""
    tree = etree.parse(os.path.join(os.path.dirname(os.path.realpath(__file__)), "template.xml"))

    attstrCategory = tree.xpath('//params/section[@name="Header"]/attstr[@name="category"]')[0]
    attstrCategory.set("val", parameters["category"])

    attstrAuthor = tree.xpath('//params/section[@name="Header"]/attstr[@name="author"]')[0]
    if parameters["author"] == "":
        attstrAuthor.getparent().remove(attstrAuthor)
    else:
        attstrAuthor.set("val", parameters["author"])

    attstrName = tree.xpath('//params/section[@name="Header"]/attstr[@name="name"]')[0]
    attstrName.set("val", parameters["name"])

    attstrDescription = tree.xpath('//params/section[@name="Header"]/attstr[@name="description"]')[0]
    attstrDescription.set("val", parameters["description"])

    attstrGraphic = tree.xpath('//params/section[@name="Graphic"]/attstr[@name="3d description"]')[0]
    attstrGraphic.set("val", os.path.splitext(os.path.basename(parameters["output_file"]))[0] + ".ac")

    attnumStepLength = tree.xpath('//params/section[@name="Main Track"]/attnum[@name="profil steps length"]')[0]
    attnumStepLength.set("val", str(parameters["step_length"]))

    node = tree.xpath('//params/section[@name="Main Track"]/section[@name="Track Segments"]')[0]

    # Remove all child elements of the found node
    for child in list(node):
        node.remove(child)

    tracks = generate_segments_xml(segments, parameters)
    parent = node.getparent()
    parent.replace(node, tracks)

    dump = etree.tostring(tree, encoding="UTF-8", xml_declaration=True)
    formatter = xmlformatter.Formatter(indent="2", selfclose=True)
    pretty_xml_as_bytes = formatter.format_string(dump)
    pretty_xml_as_string = pretty_xml_as_bytes.decode("utf-8")

    dump = pretty_xml_as_string.replace("&amp;default-surfaces;", "&default-surfaces;")
    with open(parameters["output_file"], "w") as file:
        file.write(dump)
    

def format_number(value):
    return f"{value:.10f}".rstrip("0").rstrip(".")

def generate_segments_xml(segments, parameters):
    tracks = etree.Element("section")
    tracks.set("name", "Track Segments")
    section_number = 0
    scale = parameters["track_length"] / parameters["svg_length"] if parameters["track_length"] != 0 else 1.0
    start_offset = parameters["start_offset"] / scale
    named_pieces = []
    for i, segment in enumerate(segments):
        is_curved = not isinstance(segment, Line)

        if i % 2 == 1 and not is_curved:
            raise Exception("The track must have a curve after a straight segment")

        # Each piece is (angle, radius, end radius, length); straights have no radius
        if is_curved:
            segment_name = f"Curve {(i + 1) / 2:.0f}"
            prev_segment = segments[i - 1]
            next_segment = segments[(i + 1) % len(segments)]

            angle = angle_between_segments(prev_segment, next_segment)
            if abs(abs(angle) - np.pi) < 1e-9:
                # A U-turn between parallel straights has no sign of its own: turn towards its end
                angle = np.copysign(np.pi, np.angle((segment.end - segment.start) / _heading(prev_segment)))
            radii = calculate_radius_end_radius(angle, segment, prev_segment, scale, parameters["step_length"])
            if radii is not None:
                pieces = [(angle, radii[0], radii[1], segment.length())]
            else:
                pieces = calculate_biarc(segment, prev_segment, next_segment)
                print(f"Warning: {segment_name} cannot be built with a single curve, it was replaced by: "
                      + " + ".join(describe_piece(piece, parameters) for piece in pieces))
        else:
            segment_name = f"Straight {(i + 2) / 2:.0f}"
            pieces = [(0, None, None, segment.length())]
            if i == 0 and start_offset > 0:
                # The track starts at the start line, so the first straight begins there
                pieces = [(0, None, None, segment.length() - start_offset)]

        named_pieces.append((segment_name, pieces))

    if start_offset > 0:
        # The part of the first straight that is before the start line closes the track
        named_pieces.append(("Straight 1 (before start)", [(0, None, None, start_offset)]))

    for segment_name, pieces in named_pieces:
        for n, (angle, radius, end_radius, length) in enumerate(pieces):
            if radius is None:
                curvature = "str"
            elif angle > 0:
                curvature = "rgt"
            else:
                curvature = "lft"

            piece_name = segment_name if len(pieces) == 1 else f"{segment_name} ({n + 1}/{len(pieces)})"
            tracks.append(etree.Comment(f"******************************"))
            tracks.append(etree.Comment(f"{piece_name.center(30)}"))
            tracks.append(etree.Comment(f"******************************"))

            section_number += 1
            track = etree.SubElement(tracks, "section")
            track.set("name", f'{section_number}')
            attstr = etree.SubElement(track, "attstr")
            attstr.set("name", "type")
            attstr.set("val", curvature)

            if curvature == "str":
                attnum = etree.SubElement(track, "attnum")
                attnum.set("name", "lg")
                attnum.set("unit", "m")
                attnum.set("val", length_in_meters(length, parameters))
            else:
                attnum = etree.SubElement(track, "attnum")
                attnum.set("name", "arc")
                attnum.set("unit", "deg")
                attnum.set("val", format_number(np.degrees(abs(angle))))
                attnum = etree.SubElement(track, "attnum")
                attnum.set("name", "radius")
                attnum.set("unit", "m")
                attnum.set("val", length_in_meters(radius, parameters))
                if end_radius != radius:
                    attnum = etree.SubElement(track, "attnum")
                    attnum.set("name", "end radius")
                    attnum.set("unit", "m")
                    attnum.set("val", length_in_meters(end_radius, parameters))
            attstr = etree.SubElement(track, "attstr")
            attstr.set("name", "surface")
            attstr.set("val", parameters["surface"])

    return tracks

def describe_piece(piece, parameters):
    angle, radius, end_radius, length = piece
    if radius is None:
        return f"straight of {float(length_in_meters(length, parameters)):.1f} m"
    return f"{'right' if angle > 0 else 'left'} curve of {np.degrees(abs(angle)):.2f} deg with radius {float(length_in_meters(radius, parameters)):.1f} m"

def get_output_file(parameters, parser):
    if not os.path.exists(parameters["track_directory"]):
        parser.error(f"The track directory '{parameters['track_directory']}' does not exist.")

    output_file = os.path.join(parameters["track_directory"], parameters["category"])
    if not os.path.exists(output_file):
        try:
            os.makedirs(output_file)
            print(f"Creating track directory: {output_file}")
        except Exception as e:
            parser.error(f"Failed to create track directory '{output_file}'.")

    codename = parameters["name"].lower().replace(" ", "_")
    output_file = os.path.join(output_file, codename)
    if not os.path.exists(output_file):
        try:
            os.makedirs(output_file)
            print(f"Creating track directory: {output_file}")
        except Exception as e:
            parser.error(f"Failed to create track directory '{output_file}'.")
    
    output_file = os.path.join(output_file, codename) + ".xml"

    return output_file

def main():
    # Object with the parameters
    parameters = {
        "svg_length": 0,
        "track_length": 0,
        "track_directory": "",
        "surface": "",
        "step_length": 4.0,
        "start_offset": 0.0,
        "category": "",
        "author": "",
        "name": "",
        "description": "",
    }

    if os.name == "posix":
        parameters["default_track_directory"] = "/usr/share/games/torcs/tracks"
    elif os.name == "nt":
        parameters["default_track_directory"] = "C:\\Program Files\\Torcs\\tracks"
    elif os.name == "darwin":
        parameters["default_track_directory"] = "/Applications/Torcs.app/Contents/Resources/tracks"
    else:
        raise ValueError(f"Unsupported operating system: {os.name}")

    parser = argparse.ArgumentParser(description="Convert SVG file to XML track for TrackGen for Torcs and Speed Dreams")
    parser.add_argument("svgfile", type=str, help="Path of the SVG file")
    parser.add_argument("-t", "--track-directory", type=str, default=parameters["default_track_directory"], help="The track directory of Torcs, Speed Dreams, or any other destination directory. Default: " + parameters["default_track_directory"])
    parser.add_argument("-c", "--category", type=str, default="circuit", help="The track category (circuit, road, oval, dirt...). Default is 'circuit'")
    parser.add_argument("-l", "--length", type=float, default=0, help="The length of the first segment in meters")
    parser.add_argument("-s", "--surface", type=str, default="asphalt", help="The surface of the track. Possible values are: 'asphalt-lines', 'asphalt-l-left', 'asphalt-l-right', 'asphalt-l-both', 'asphalt-pits', 'asphalt', 'dirt', 'dirt-b', 'asphalt2', 'road1', 'road1-pits', 'road1-asphalt', 'asphalt-road1', 'b-road1', 'b-road1-l2', 'b-road1-l2p', 'concrete', 'concrete2', 'concrete3', 'b-asphalt', 'b-asphalt-l1', 'b-asphalt-l1p', 'asphalt2-lines', 'asphalt2-l-left', 'asphalt2-l-right', 'asphalt2-l-both', 'grass', 'grass3', 'grass5', 'grass6', 'grass7', 'gravel', 'sand3', 'sand', 'curb-5cm-r', 'curb-5cm-l', 'curb-l', 'tar-grass3-l', 'tar-grass3-r', 'tar-sand'. Default is 'asphalt'")
    parser.add_argument("--step-length", type=float, default=4.0, help="The profil steps length of the track in meters. Default is 4.0")
    parser.add_argument("--start-offset", type=float, default=0, help="Distance in meters from the beginning of the first segment to the start line. Default is 0")
    parser.add_argument("-a", "--author", type=str, default="", help="Name of author of the track")
    parser.add_argument("-n", "--name", type=str, default="Oval", help="The track name")
    parser.add_argument("-d", "--description", type=str, help="The description of the track")

    
    args = parser.parse_args()

    parameters["track_directory"] = args.track_directory
    parameters["main_straight_length"] = args.length
    parameters["surface"] = args.surface
    parameters["step_length"] = args.step_length
    parameters["start_offset"] = args.start_offset
    parameters["category"] = args.category
    parameters["author"] = args.author
    parameters["name"] = args.name
    parameters["description"] = args.description

    if parameters["step_length"] <= 0:
        parser.error("The step length must be greater than 0.")

    parameters["output_file"] = get_output_file(parameters, parser)
    print(f"Output file: {parameters['output_file']}")

    # Load SVG file and get the segments of the path
    try:
        paths, attributes = svg2paths(args.svgfile)
    except Exception as e:
        parser.error(f"Failed to parse the SVG file '{args.svgfile}': {e}")
    segments = paths[0]

    # Chechk if the path is closed
    is_closed = segments.isclosed()
    print(f"The path {'is circular' if is_closed else 'is not circular'}.")

    # Rounding errors can make the closing "z" of the SVG path add a line of almost zero
    # length, which would become an empty last segment with a meaningless direction
    if is_closed and isinstance(segments[-1], Line) and segments[-1].length() < 1e-9 * segments.length():
        segments = Path(*segments[:-1])

    # Calculate the total length of the track
    parameters["svg_length"] = segments.length()
    parameters["track_length"] = parameters["main_straight_length"] / segments[0].length() * parameters["svg_length"]
    print(f"The total length of the track is " + length_in_meters(parameters['svg_length'], parameters) + " meters.")
    print(f"Total segments: {len(segments)}")

    first_segment_length = float(length_in_meters(segments[0].length(), parameters))
    if not 0 <= parameters["start_offset"] < first_segment_length:
        parser.error(f"The start offset must be between 0 and the length of the first segment ({first_segment_length:.2f} m).")

    create_xml_file(segments, parameters)

if __name__ == "__main__":
    main()
