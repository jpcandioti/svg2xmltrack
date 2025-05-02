SVG2XMLTrack
============


## Installation

### Dependencies

~~~.shell
pip install argparse svgpathtools numpy lxml xmlformatter
~~~

### Download

~~~.shell
git clone https://github.com/jpcandioti/svg2xmltrack.git
~~~

## How to create SVG file

You can create an SVG file with many applications and on different platforms (web, desktop, mobile). A good option is [Inkscape](https://inkscape.org/).
If you want to create a real circuit, you can paste a map image and draw the SVG over it.

### Requirements
- Always start with a straight segment.
- You must intersperse straight and curved segments.
- If the path is circular, the last segment must be curved.


## Usage and Options

~~~.shell
$ ./svg2xmltrack.py -h
usage: svg2xmltrack.py [-h] [-t TRACK_DIRECTORY] [-c CATEGORY] [-l LENGTH] [-s SURFACE] [-a AUTHOR] [-n NAME] [-d DESCRIPTION] svgfile

Convert SVG file to XML track for TrackGen for Torcs and Speed Dreams

positional arguments:
  svgfile               Pathame of the SVG file

optional arguments:
  -h, --help            show this help message and exit
  -t TRACK_DIRECTORY, --track-directory TRACK_DIRECTORY
                        The track directory of Torcs, Speed Dreams, or any other destination directory. Default: /usr/share/games/torcs/tracks
  -c CATEGORY, --category CATEGORY
                        The track category (circuit, road, oval, dirt...). Default is 'circuit'
  -l LENGTH, --length LENGTH
                        The length of the first segment in meters
  -s SURFACE, --surface SURFACE
                        The surface of the track. Possible values are: 'asphalt-lines', 'asphalt-l-left', 'asphalt-l-right', 'asphalt-l-both', 'asphalt-pits', 'asphalt', 'dirt', 'dirt-b', 'asphalt2', 'road1', 'road1-pits', 'road1-asphalt', 'asphalt-road1',
                        'b-road1', 'b-road1-l2', 'b-road1-l2p', 'concrete', 'concrete2', 'concrete3', 'b-asphalt', 'b-asphalt-l1', 'b-asphalt-l1p', 'asphalt2-lines', 'asphalt2-l-left', 'asphalt2-l-right', 'asphalt2-l-both', 'grass', 'grass3', 'grass5',
                        'grass6', 'grass7', 'gravel', 'sand3', 'sand', 'curb-5cm-r', 'curb-5cm-l', 'curb-l', 'tar-grass3-l', 'tar-grass3-r', 'tar-sand'. Default is 'asphalt'
  -a AUTHOR, --author AUTHOR
                        Name of author of the track
  -n NAME, --name NAME  The track name
  -d DESCRIPTION, --description DESCRIPTION
                        The description of the track
~~~


## Example

~~~.shell
./svg2xmltrack.py -nRosario -ccircuit -d"Autódromo Juan Manuel Fangio de la Ciudad de Rosario" -l480.53 rosario.svg
~~~
