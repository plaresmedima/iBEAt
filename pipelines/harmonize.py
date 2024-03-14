import numpy as np
import dbdicom

# Note from older code:
# For T1-mapping, Siemens uses the field 'InversionTime' 
# but Philips uses (0x2005, 0x1572). Need to include a 
# harmonization step for Philips.

def magnetization_transfer(database):
    # Merge MT ON and MT OFF sequences and remove the originals

    # Find the series
    mt_on = database.series(SeriesDescription="MT_ON_kidneys_cor-oblique_bh")
    mt_off = database.series(SeriesDescription="MT_OFF_kidneys_cor-oblique_bh")

    # If one is missing, or there are multiple, do not proceed
    if len(mt_on)!=1 or len(mt_off)!=1:
        return
    
    # Merge the two sequences into one
    try:
        study = mt_on[0].parent()
        mt = study.new_series(SeriesDescription='MT_kidneys_cor-oblique_bh')
        dbdicom.merge(mt_on+mt_off, mt)
    except:
        return
    
    # If successfully merged, remove the originals
    mt_on[0].remove()
    mt_off[0].remove()



def ivim(database):
    desc = 'IVIM_kidneys_cor-oblique_fb'
    series = database.series(SeriesDescription=desc)[0]

    if series.instance().Manufacturer=='SIEMENS':
        # The b-values and bvectors are incorrectly stored in the header, so must be hardwired.
        # First create an array of b-values and b-vectors ordered by slice location and acquisition time.
        # for each slice, 10 b-values are acquired 3 times, once in each direction (-x, -y, z).
        dims = ('SliceLocation', 'InstanceNumber')
        nslices = 30
        nbvals = 10
        bvals = [0, 10, 20, 30, 50, 80, 100, 200, 300, 600.0]
        bvals = np.array(nslices*[bvals + bvals + bvals])
        bvecs_z = [[-1,0,0]]*nbvals + [[0,-1,0]]*nbvals + [[0,0,1]]*nbvals
        bvecs = np.zeros((nslices,3*nbvals), dtype=object)
        for z in range(nslices):
            for b in range(3*nbvals):
                bvecs[z,b] = bvecs_z[b]
        # Write b-values and b-vectors in the correct DICOM data element
        series.set_values((bvals,bvecs), ('DiffusionBValue', 'DiffusionGradientOrientation'), dims=dims)

def dti(database):
    desc = "DTI_kidneys_cor-oblique_fb"
    series = database.series(SeriesDescription=desc)[0]

    if series.instance().Manufacturer=='SIEMENS':
        # bvalues and b-vectors are correctly stored but in a private data tag
        # Copy the values to the standard DICOM data tag
        bvals, bvecs = series.values((0x19, 0x100c), (0x19, 0x100e)) 
        series.set_values((bvals,bvecs), ('DiffusionBValue', 'DiffusionGradientOrientation'))


def t2(database):
    desc = "T2map_kidneys_cor-oblique_mbh_magnitude"
    series = database.series(SeriesDescription=desc)[0]

    if series.instance().Manufacturer=='SIEMENS':
        TE = [[0,30,40,50,60,70,80,90,100,110,120]]*5
        series.set_values(np.array(TE), 'InversionTime',  dims=('SliceLocation', 'InstanceNumber'))


def dce(database):
    desc = "DCE_kidneys_cor-oblique_fb"
    series = database.series(SeriesDescription=desc)[0]

    if series.instance().Manufacturer=='SIEMENS':
        # Split into two separate series, one for the aorta and one for the kidneys
        # Then delete the original
        try:
            split_series = series.split_by('ImageOrientationPatient')
        except:
            # If the series has already been split, there is nothing to be done.
            return
        locs = split_series[0].unique('SliceLocation')
        if len(locs) == 1:
            split_series[0].SeriesDescription = "DCE_aorta_axial_fb"
            split_series[1].SeriesDescription = desc
        else:
            split_series[0].SeriesDescription = desc
            split_series[1].SeriesDescription = "DCE_aorta_axial_fb"
        series.remove()
        	



def all_series(database):

    # magnetization_transfer(database)
    # ivim(database)
    # dti(database)
    # dce(database)
    t2(database)


