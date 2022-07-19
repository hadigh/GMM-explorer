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
import openquake.hazardlib.imt as hi
import numpy as np
import pandas as pd

# inputs
db_file = '/home/hadi/gmprocess_projects/data/default_default_metrics_rotd50.0.csv'

# output
output_file = '../results/db_with_gmms.csv'
# read db
db = pd.read_csv(db_file)
db = db[(db['HypocentralDistance'] <= 300) & (db['EarthquakeMagnitude'] >= 4.0)]
db.reset_index(drop=True, inplace=True)
# db = db.head(30)
# get list of GMMs
gs = get_available_gsims()
# gmms = [gs['Allen2012'](), gs['AbrahamsonSilva2008'](), gs['ChiouYoungs2008'](), gs['CampbellBozorgnia2008']()] #avoids importing classes, G = 3
# gmms = [gs['Allen2012']()]

# for gmm in gmms:
# for n , gg in get_available_gsims().items():
#     ...:     aa = gg.DEFINED_FOR_TECTONIC_REGION_TYPE
#     ...:     try:
#     ...:         print(aa.value)
#     ...:     except AttributeError:
#     ...:         print(aa)
passed_gmm = []
failed_gmm = []
for name, gmm in get_available_gsims().items():

    trt = gmm.DEFINED_FOR_TECTONIC_REGION_TYPE
    imt = gmm.DEFINED_FOR_INTENSITY_MEASURE_TYPES

    try:

        if ((trt.name == 'ACTIVE_SHALLOW_CRUST') | (trt.name == 'STABLE_CONTINENTAL')) & (hi.SA in imt) & \
                (not ((name == 'ChiouYoungs2014ACME2019') | (name == 'ChiouYoungs2014NearFaultEffect')
                      | (name == 'GulerceAbrahamson2011') | (name == 'MegawatiEtAl2003')
                      | (name == 'SgobbaEtAl2020') | (name == 'YuEtAl2013Ms')
                      | (name == 'YuEtAl2013MsEastern') | (name == 'YuEtAl2013MsStable')
                      | (name == 'YuEtAl2013MsTibet') | (name == 'YuEtAl2013Mw')
                      | (name == 'YuEtAl2013MwEastern') | (name == 'YuEtAl2013MwStable')
                      | (name == 'YuEtAl2013MwTibet')

                )):
            
            ####(OQ3.13: | (name =='AmeriEtAl2017RepiStressDrop')
                      # | (name == 'Bradley2013bChchMaps') | (name == 'Bradley2013bChchMapsAdditionalSigma')

            db['_'.join((name, 'mean', 'SA(0.1)'))] = np.nan
            db['_'.join((name, 'std', 'SA(0.1)'))] = np.nan

            db['_'.join((name, 'mean', 'SA(1)'))] = np.nan
            db['_'.join((name, 'std', 'SA(1)'))] = np.nan

            # dum_med_sa

            for index, row in db.iterrows():

                # this is the site for which we compute the median ground shaking
                sta_loc = Point(db['StationLongitude'][index], db['StationLatitude'][index])
                site_collection = SiteCollection([Site(location=sta_loc, vs30=760., vs30measured=False, z1pt0=40., z2pt5=1.0)])

                # this is the source for which we compute the median ground shaking
                src_loc = Point(db['EarthquakeLongitude'][index], db['EarthquakeLatitude'][index])
                src_depth = db['EarthquakeDepth'][index]

                src_depth = 0.0 if src_depth < 0.0 else src_depth # some of the events have negative depth! check with Trev (2001-10-09T15:56:05.000000Z)


                src_mag = db['EarthquakeMagnitude'][index]
                src = PointSource(
                    source_id=str(db['EarthquakeId'][index]),
                    name='point',
                    tectonic_region_type='Active Shallow Crust',
                    mfd=TruncatedGRMFD(min_mag=src_mag, max_mag=src_mag + 0.1, bin_width=0.1, a_val=0.01, b_val=0.98),
                    rupture_mesh_spacing=2.,
                    magnitude_scaling_relationship=WC1994(),
                    rupture_aspect_ratio=1.,
                    temporal_occurrence_model=PoissonTOM(50.),
                    upper_seismogenic_depth=0.,
                    lower_seismogenic_depth=50.,
                    location=src_loc,
                    nodal_plane_distribution=PMF([(1., NodalPlane(strike=45, dip=50, rake=0))]),
                    hypocenter_distribution=PMF([(1, src_depth)]) #depth in km
                )
                ruptures = [r for r in src.iter_ruptures()] #create list of ruptures
                ruptures[0].mag = src_mag # this is to make sure the magnitude is consistent and avoid python precision!

                # these are the intensity measure type for which we compute the median ground shaking
                imtls = {s: [0] for s in ['SA(0.1)', 'SA(1)']}  # required for context maker, M = 2 IMTs

                context_maker = ContextMaker('*', [gmm()], {'imtls': imtls})  # necessary contexts builder
                ctxs = context_maker.get_ctxs(ruptures, site_collection)  # returns rupture contexts
                gms = context_maker.get_mean_stds(
                    ctxs)  # calculate ground motions and stds, returns array of shape (4, G, M, N)

                db.loc[index, '_'.join((name, 'mean', 'SA(0.1)'))] = gms[0][0][0][0]
                db.loc[index, '_'.join((name, 'std', 'SA(0.1)'))] = gms[1][0][0][0]

                db.loc[index, '_'.join((name, 'mean', 'SA(1)'))] = gms[0][0][1][0]
                db.loc[index, '_'.join((name, 'std', 'SA(1)'))] = gms[1][0][1][0]
                # import pdb; pdb.set_trace()
                # print(gmm())
                # print(gms)

                # if (gms.size==0):
                if np.isnan(gms[0][0][1][0]):
                    print("OH No: %s %f" % (name, db['EpicentralDistance'][index]))
                    import pdb; pdb.set_trace()

            passed_gmm.append(name)

    except AttributeError:
        print("Could not %s" %name)
        failed_gmm.append(name)

        pass

db.to_csv(output_file)

failed_gmm = {'failed_gmms': failed_gmm}
failed_gmm = pd.DataFrame(failed_gmm)
failed_gmm.to_csv("../results/failed_gmms.csv")

passed_gmm = {'passed_gmms': passed_gmm}
passed_gmm = pd.DataFrame(passed_gmm)
passed_gmm.to_csv("../results/passed_gmm.csv")

#

#
# # these are the intensity measure type for which we compute the median ground shaking
# imtls = {s: [0] for s in ['PGA','SA(0.3)']} #required for context maker, M = 2 IMTs
#
# context_maker = ContextMaker('*',gmpes,{'imtls': imtls}) #necessary contexts builder
# ctxs = context_maker.get_ctxs(ruptures,site_collection) #returns rupture contexts
# gms = context_maker.get_mean_stds(ctxs) #calculate ground motions and stds, returns array of shape (4, G, M, N)
# print(gms.shape) #first 4 are median, std_total, std_intra, std_inter, then G=3 gsims, M=2 IMTs, 15 scenarios = magnitudes
#
# print(np.exp(gms[0][2][0])) #median values, AbrahamsonSilva2008, PGA