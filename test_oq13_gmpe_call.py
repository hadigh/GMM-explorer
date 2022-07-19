#to switch to OQ3.13 (run the following once in ipython
# import sys
# sys.path.append('/home/hadi/openquake/lib/python3.9/site-packages/openquake/hazardlib/')


import os
import pandas as pd
# from openquake.hazardlib.site import Site, SiteCollection
from geo import Point
from sourceconverter import SourceConverter
from nrml import to_python
# from gsim.abrahamson_2014 import AbrahamsonEtAl2014
# from gsim.allen_2012 import Allen2012
# # from openquake.hazardlib.gsim.chiou_youngs_2014 import ChiouYoungs2014ACME2019
# from gsim.kotha_2020 import KothaEtAl2020
from gsim.nga_east import NGAEastGMPE
from contexts import ContextMaker
from imt import PGA, SA
from const import StdDev
from source import PointSource
from mfd import TruncatedGRMFD
from scalerel import WC1994
from geo import Point, NodalPlane, Line
from pmf import PMF
from tom import PoissonTOM
from gsim.gmpe_table import GMPETable



fname = os.path.join('../data/sites.csv')
sites = pd.read_csv(fname, names=['lon', 'lat'])

# Create site-collection
sites_list = []
for i, row in sites.iterrows():
    site = Site(Point(row.lon, row.lat), vs30=760, z1pt0=30, z2pt5=0.5,
                vs30measured=True)
    sites_list.append(site)
sitec = SiteCollection(sites_list)

# Create source
location = Point(-120., 35.)
src = PointSource(
    source_id='1',
    name='point',
    tectonic_region_type='Active Shallow Crust',
    mfd=TruncatedGRMFD(min_mag=5., max_mag=6.5, bin_width=0.1, a_val=0.01, b_val=0.98),
    rupture_mesh_spacing=2.,
    magnitude_scaling_relationship=WC1994(),
    rupture_aspect_ratio=1.,
    temporal_occurrence_model=PoissonTOM(50.),
    upper_seismogenic_depth=2.,
    lower_seismogenic_depth=12.,
    location=location,
    nodal_plane_distribution=PMF([(1., NodalPlane(strike=45, dip=50, rake=0))]),
    hypocenter_distribution=PMF([(1, 7.)])
)
rup = [r for r in src.iter_ruptures()][0]
# Create the contexts
gmm = AbrahamsonEtAl2014()
# gmm = Allen2012()
# gmm = ChiouYoungs2014ACME2019()
# gmm = KothaEtAl2020()
# gmm = NGAEastGMPE(gmpe_table='NGAEast_YENIER_ATKINSON.hdf5')
param = dict(imtls={'PGA': []}, cache_distances=True)
cm = ContextMaker('*', [gmm], param)
[ctx] = cm.get_ctxs([rup], sitec)
imt = PGA()
gmm_mean, gmm_std, = gmm.get_mean_and_stddevs(ctx,ctx, ctx, imt, [])

print(gmm_mean, gmm_std)
