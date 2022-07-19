from openquake.hazardlib.gsim import get_available_gsims #building this takes time
from openquake.hazardlib.source import PointSource
from openquake.hazardlib.mfd import TruncatedGRMFD
from openquake.hazardlib.scalerel import WC1994
from openquake.hazardlib.geo import Point, NodalPlane
from openquake.hazardlib.pmf import PMF
from openquake.hazardlib.tom import PoissonTOM
from openquake.hazardlib.site import Site, SiteCollection
from openquake.hazardlib.imt import PGA
from openquake.hazardlib.const import StdDev
from openquake.hazardlib.contexts import ContextMaker
import numpy as np


gs = get_available_gsims()
gmpes = [gs['AbrahamsonSilva2008'](), gs['ChiouYoungs2008'](), gs['CampbellBozorgnia2008'](),] #avoids importing classes, G = 3


# explore magnitude scaling, by defining a Point source and calculating median ground shaking at the point source location
location = Point(9.1500, 45.1833)
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
    hypocenter_distribution=PMF([(1, 10.)]) #depth 10 km
)
ruptures = [r for r in src.iter_ruptures()] #create list of ruptures
mags = [r.mag for r in ruptures] #magnitude list - from 5.05 to 6.45 step 0.1 (15)

# this is the site for which we compute the median ground shaking
site_collection = SiteCollection([Site(location=location, vs30=760., vs30measured=True, z1pt0=40., z2pt5=1.0)])

# these are the intensity measure type for which we compute the median ground shaking
imtls = {s: [0] for s in ['PGA','SA(0.3)']} #required for context maker, M = 2 IMTs

context_maker = ContextMaker('*',gmpes,{'imtls': imtls}) #necessary contexts builder
ctxs = context_maker.get_ctxs(ruptures,site_collection) #returns rupture contexts
gms = context_maker.get_mean_stds(ctxs) #calculate ground motions and stds, returns array of shape (4, G, M, N)
print(gms.shape) #first 4 are median, std_total, std_intra, std_inter, then G=3 gsims, M=2 IMTs, 15 scenarios = magnitudes

print(np.exp(gms[0][2][0])) #median values, AbrahamsonSilva2008, PGA