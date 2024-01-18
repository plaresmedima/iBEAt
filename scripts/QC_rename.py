import os
import pandas as pd

def main(database):
    df = pd.DataFrame([
        ['T2w_abdomen_haste_tra_mbh',0],
        ['T1w_abdomen_dixon_cor_bh_out_phase',0],
        ['T1w_abdomen_dixon_cor_bh_in_phase',0],
        ['T1w_abdomen_dixon_cor_bh_fat',0],
        ['T1w_abdomen_dixon_cor_bh_water',0],
        ['PC_RenalArtery_Right_EcgTrig_fb_120',0],
        ['PC_RenalArtery_Right_EcgTrig_fb_120_magnitude',0],
        ['PC_RenalArtery_Right_EcgTrig_fb_120_phase',0],
        ['PC_RenalArtery_Left_EcgTrig_fb_120',0],
        ['PC_RenalArtery_Left_EcgTrig_fb_120_magnitude',0],
        ['PC_RenalArtery_Left_EcgTrig_fb_120_phase',0],
        ['T2star_map_pancreas_tra_mbh_magnitude',0],
        ['T2star_map_pancreas_tra_mbh_phase',0],
        ['T2star_map_pancreas_tra_mbh_T2star',0],
        ['T1w_kidneys_cor-oblique_mbh_magnitude',0],
        ['T1w_kidneys_cor-oblique_mbh_phase',0],
        ['T1map_kidneys_cor-oblique_mbh_magnitude',0],
        ['T1map_kidneys_cor-oblique_mbh_phase',0],
        ['T1map_kidneys_cor-oblique_mbh_moco',0],
        ['T1map_kidneys_cor-oblique_mbh_T1map',0],
        ['T2map_kidneys_cor-oblique_mbh_magnitude',0],
        ['T2map_kidneys_cor-oblique_mbh_phase',0],
        ['T2map_kidneys_cor-oblique_mbh_moco',0],
        ['T2map_kidneys_cor-oblique_mbh_T2map',0],
        ['T2star_map_kidneys_cor-oblique_mbh_magnitude',0],
        ['T2star_map_kidneys_cor-oblique_mbh_phase',0],
        ['T2star_map_kidneys_cor-oblique_mbh_T2star',0],
        ['DTI_kidneys_cor-oblique_fb',0],
        ['IVIM_kidneys_cor-oblique_fb',0],
        ['MT_OFF_kidneys_cor-oblique_bh',0],
        ['MT_ON_kidneys_cor-oblique_bh',0],
        ['ASL_kidneys_pCASL_cor-oblique_fb_RBF_moco',0],
        ['DCE_kidneys_cor-oblique_fb',0],
        ['T1w_abdomen_dixon_cor_bh_out_phase_post_contrast',0],
        ['T1w_abdomen_dixon_cor_bh_in_phase_post_contrast',0],
        ['T1w_abdomen_dixon_cor_bh_fat_post_contrast',0],
        ['T1w_abdomen_dixon_cor_bh_water_post_contrast',0]
        ], columns=['MRI Sequence','Checked'])

    list_of_series = database.series()
    list_of_series_description = []
    for series in (list_of_series):
        list_of_series_description.append(series['SeriesDescription'])
    
    # Loop through the rows and update the second column based on the condition
    for index, row in df.iterrows():
        if row['MRI Sequence'] in list_of_series_description:
            df.at[index, 'Checked'] = 1

    def color_rule(val):
        return ['background-color: red' if x == 0 else 'background-color: green' for x in val]

    iBEAt_column = df.style.apply(color_rule, axis=1, subset=['Checked'])
    iBEAt_column.to_excel(os.path.join(database.path(),'QC','iBEAt_checklist.xlsx'), engine='openpyxl', index=False)

    



