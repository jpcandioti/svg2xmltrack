#!/usr/bin/env python3

import argparse
import os
from svgpathtools import svg2paths
from svgpathtools.path import Line
import numpy as np
from lxml import etree
import xmlformatter


def calculate_radius_end_radius(angle_radians, segment, prev_segment, next_segment):
    angle_radians = abs(angle_radians)
    distance = abs(next_segment.start - prev_segment.end)

    angle_real = abs(angle_between_segments(prev_segment, segment) * 2)

    distorsion = angle_real - angle_radians

    angle_end = abs(angle_radians - distorsion)

    distance_start = distance / angle_radians * angle_end
    distance_end = distance / angle_radians * angle_real

    radius_start = distance_start / (2 * np.sin(angle_radians / 2))
    radius_end = distance_end / (2 * np.sin(angle_radians / 2))
   
    return radius_start, radius_end

def angle_between_segments(segment1, segment2):
    vector1 = np.array([segment1.end.real - segment1.start.real, segment1.end.imag - segment1.start.imag])
    vector2 = np.array([segment2.end.real - segment2.start.real, segment2.end.imag - segment2.start.imag])

    unit_v1 = vector1 / np.linalg.norm(vector1)
    unit_v2 = vector2 / np.linalg.norm(vector2)
    dot_product = np.dot(unit_v1, unit_v2)

    angle_radians = np.arccos(np.clip(dot_product, -1.0, 1.0))

    cross_product = np.cross(vector1, vector2)

    if cross_product < 0:
        angle_radians = -angle_radians

    return angle_radians

def length_in_meters(value, parameters):
    if parameters["track_length"] == 0:
        return value

    return format_number((value / parameters["svg_length"]) * parameters["track_length"])

def create_xml_file(segments, parameters):
    doctype = """<!DOCTYPE params SYSTEM "../../../src/libs/tgf/params.dtd" [
<!--  general definitions for tracks  -->
<!ENTITY default-surfaces SYSTEM "../../../data/tracks/surfaces.xml">
]>"""
    tree = etree.parse("template.xml")

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
    attstrGraphic.set("val", parameters["name"] + ".ac")

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
    for i, segment in enumerate(segments):
        length = segment.length()
        is_curved = not isinstance(segment, Line)

        if is_curved:
            prev_segment = segments[i - 1]
            next_segment = segments[(i + 1) % len(segments)]

            angle = angle_between_segments(prev_segment, next_segment)
            if angle > 0:
                curvature = "rgt"
            else:
                curvature = "lft"

            radius, end_radius = calculate_radius_end_radius(angle, segment, prev_segment, next_segment)
        else:
            curvature = "str"

        if i % 2 == 1 and curvature == "str":
            raise Exception("The track must have a curve after a straight segment")

        if is_curved:
            segment_name = f"Curve {(i + 1) / 2:.0f}"
        else:
            segment_name = f"Straight {(i + 2) / 2:.0f}"

        tracks.append(etree.Comment(f"******************************"))
        tracks.append(etree.Comment(f"     {segment_name}                "))
        tracks.append(etree.Comment(f"******************************"))

        track = etree.SubElement(tracks, "section")
        track.set("name", f'{i + 1}')
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

def get_output_file(parameters, parser):
    if not os.path.exists(parameters["track_directory"]):
        parser.error(f"The track directory '{parameters['track_directory']}' does not exist.")

    output_file = os.path.join(parameters["track_directory"], parameters["category"])
    if not os.path.exists(output_file):
        try:
            os.makedirs(output_file)
            print(f"Creating track directory: {parameters['output_file']}")
        except Exception as e:
            parser.error(f"Failed to create track directory '{parameters['output_file']}'.")

    codename = parameters["name"].lower().replace(" ", "_")
    output_file = os.path.join(output_file, codename)
    if not os.path.exists(output_file):
        try:
            os.makedirs(output_file)
            print(f"Creating track directory: {parameters['output_file']}")
        except Exception as e:
            parser.error(f"Failed to create track directory '{parameters['output_file']}'.")
    
    output_file = os.path.join(output_file, codename) + ".xml"

    return output_file

def main():
    # Object with the parameters
    parameters = {
        "svg_length": 0,
        "track_length": 0,
        "track_directory": "",
        "surface": "",
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

    parser = argparse.ArgumentParser(description="Convert SVG file to XML track for TrackGen or TrackEditor")
    parser.add_argument("svgfile", type=str, help="Pathame of the SVG file")
    parser.add_argument("-t", "--track-directory", type=str, default=parameters["default_track_directory"], help="The track directory of Torcs, Speed Dreams, or any other destination directory. Default: " + parameters["default_track_directory"])
    parser.add_argument("-c", "--category", type=str, default="circuit", help="The track category (circuit, road, oval, dirt...). Default is 'circuit'")
    parser.add_argument("-l", "--length", type=float, default=0, help="The length of the first segment in meters")
    parser.add_argument("-s", "--surface", type=str, default="asphalt", help="The surface of the track. Possible values are: 'asphalt-lines', 'asphalt-l-left', 'asphalt-l-right', 'asphalt-l-both', 'asphalt-pits', 'asphalt', 'dirt', 'dirt-b', 'asphalt2', 'road1', 'road1-pits', 'road1-asphalt', 'asphalt-road1', 'b-road1', 'b-road1-l2', 'b-road1-l2p', 'concrete', 'concrete2', 'concrete3', 'b-asphalt', 'b-asphalt-l1', 'b-asphalt-l1p', 'asphalt2-lines', 'asphalt2-l-left', 'asphalt2-l-right', 'asphalt2-l-both', 'grass', 'grass3', 'grass5', 'grass6', 'grass7', 'gravel', 'sand3', 'sand', 'curb-5cm-r', 'curb-5cm-l', 'curb-l', 'tar-grass3-l', 'tar-grass3-r', 'tar-sand'. Default is 'asphalt'")
    parser.add_argument("-a", "--author", type=str, default="", help="Name of author of the track")
    parser.add_argument("-n", "--name", type=str, default="Oval", help="The track name")
    parser.add_argument("-d", "--description", type=str, help="The description of the track")

    
    args = parser.parse_args()

    parameters["track_directory"] = args.track_directory
    parameters["main_straight_length"] = args.length
    parameters["surface"] = args.surface
    parameters["category"] = args.category
    parameters["author"] = args.author
    parameters["name"] = args.name
    parameters["description"] = args.description

    parameters["output_file"] = get_output_file(parameters, parser)
    print(f"Output file: {parameters['output_file']}")

    # Load SVG file and get the segments of the path
    paths, attributes = svg2paths(args.svgfile)
    segments = paths[0]

    # Chechk if the path is closed
    is_closed = segments.isclosed()
    print(f"The path {'is circular' if is_closed else 'is not circular'}.")

    # Calculate the total length of the track
    parameters["svg_length"] = segments.length()
    parameters["track_length"] = parameters["main_straight_length"] / segments[0].length() * parameters["svg_length"]
    print(f"The total length of the track is " + length_in_meters(parameters['svg_length'], parameters) + " meters.")
    print(f"Total segments: {len(segments)}")

    create_xml_file(segments, parameters)

if __name__ == "__main__":
    main()
