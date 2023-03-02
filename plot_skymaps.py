# Copyright (C) 2018 Leo Singer, Michael Pürrer
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

# Please follow the installation instructions on
# https://pypi.org/project/ligo.skymap/.

from astropy.coordinates import SkyCoord
from ligo.skymap.io import fits
from ligo.skymap import plot
from ligo.skymap import postprocess
from astropy.time import Time
from astropy import units as u
import ligo.skymap.plot
import os
import json
from collections import defaultdict

import matplotlib
matplotlib.use('agg')
import matplotlib as mpl
from matplotlib import pyplot as plt


"""
Plot skymaps for O1 & O2 catalog events.
"""

def make_skymaps(fits_dict, event_colors, filename, gpstime,
                 contour_lvl=(50, 90), dpi=300, geo=True, label=True,
                 label_dicts=[dict(), dict()]):
    fig = plt.figure(figsize=(4, 4), dpi=dpi)
    # Set time
    ax = plt.axes(projection=('geo degrees mollweide' if geo
                              else 'astro hours mollweide'),
                  obstime=Time(gpstime, format='gps').utc.isot)

    # Plot contours
    for event, fits_file in fits_dict.items():
        color = event_colors[event]
        skymap, metadata = fits.read_sky_map(fits_file, nest=None)
        cls = 100 * postprocess.find_greedy_credible_levels(skymap)
        cs = ax.contour_hpx(
            (cls, 'ICRS'), nested=metadata['nest'],
            colors=color, linewidths=0.5, levels=contour_lvl)

        if label:
            # Set labels manually
            manual_labels, label_rotations = label_dicts
            for loc, rot in zip(manual_labels[event], label_rotations[event]):
                plt.text(loc[0], loc[1], event, fontsize=6, color=color, rotation=rot)

    ax.grid()

    if geo:
        # Add continents.
        geojson_filename = os.path.join(
            os.path.dirname(plot.__file__), 'ne_simplified_coastline.json')
        with open(geojson_filename, 'r') as geojson_file:
            geoms = json.load(geojson_file)['geometries']
        verts = [coord for geom in geoms
            for coord in zip(*geom['coordinates'])]
        plt.plot(*verts, color='0.5', linewidth=0.5,
            transform=ax.get_transform('world'))


    plt.savefig(filename+'.pdf')
    plt.savefig(filename+'.png')
    plt.close()


#------------------------------------------------------------------------------
if __name__ == "__main__":
    fits_dict = defaultdict(str,
                {'GW150914': 'GW150914_skymap.fits.gz',
                 'GW151226': 'GW151226_skymap.fits.gz',
                 'GW170104': 'GW170104_skymap.fits.gz',
                 'GW170608': 'GW170608_skymap.fits.gz',
                 'GW170809': 'GW170809_skymap.fits.gz',
                 'GW170814': 'GW170814_skymap.fits.gz',
                 'GW170817': 'GW170817_skymap.fits.gz',
                 'GW170818': 'GW170818_skymap.fits.gz',
                 'GW170823': 'GW170823_skymap.fits.gz',
                 'GW151012': 'GW151012_skymap.fits.gz',
                 'GW170729': 'GW170729_skymap.fits.gz'})

    event_colors = defaultdict(str,
                {'GW150914': '#00a7f0',
                 'GW151226': '#c59700',
                 'GW170104': '#55b300',
                 'GW170608': '#d48d6f',
                 'GW170809': '#00b1a4',
                 'GW170814': '#58ae87',
                 'GW170817': '#ff6c91',
                 'GW170818': '#5ca9b8',
                 'GW170823': '#9e94e2',
                 'GW151012': '#9da355',
                 'GW170729': '#ea65ff'})

    # These labels are for dpi=600
    manual_labels = { 
      'GW170104': [(40, 600), (850, 830), (1220, 395)],
      'GW170823': [(605, 345), (1500, 750)],
      'GW170608': [(1180, 750)],
      'GW151012': [(635, 370), (1310, 710)],
      'GW150914': [(1103, 150)],
      'GW170729': [(320, 400), (950, 680)],
      'GW151226': [(780, 300), (1550, 700)],
      'GW170818': [(50, 600)],
      'GW170817': [(800, 330)],
      'GW170809': [(1600, 475)],
      'GW170814': [(1150, 110)]
    }

    label_rotations = {
      'GW170104': [70, -40, -70],
      'GW170823': [70, -45],
      'GW170608': [-55],
      'GW151012': [70, -45],
      'GW150914': [15],
      'GW170729': [85, -42],
      'GW151226': [70, -38],
      'GW170818': [35],
      'GW170817': [45],
      'GW170809': [50],
      'GW170814': [12]
    }

    mpl.rcParams['font.size'] = 6

    # Confidently detected O2 GW events for which alerts were sent to EM observers
    events_in_O2_em_followup = ['GW170817', 'GW170104', 'GW170823', 'GW170608', 'GW170809', 'GW170814']
    events_all = set(fits_dict.keys())
    # O1 events (GW150914, GW151226, GW151012), along with O2 events 
    # (GW170729, GW170818) not previously released to EM observers.
    events_remaining = events_all - set(events_in_O2_em_followup)

    fits_dict1 = {k:fits_dict[k] for k in events_in_O2_em_followup}
    fits_dict2 = {k:fits_dict[k] for k in events_remaining}

    # Set time
    gpstime = 1126259462.431 # GW150914 gps time

    # plot events split into two sets in two plots
    make_skymaps(fits_dict1, event_colors, 'skymaps_O1_O2_astro_label_set1', gpstime,
                     contour_lvl=(50.0, 90.0), dpi=600, geo=False, label=True,
                     label_dicts=[manual_labels, label_rotations])
    make_skymaps(fits_dict2, event_colors, 'skymaps_O1_O2_astro_label_set2', gpstime,
                     contour_lvl=(50.0, 90.0), dpi=600, geo=False, label=True, 
                     label_dicts=[manual_labels, label_rotations])