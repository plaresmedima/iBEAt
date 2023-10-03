import os
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from io import BytesIO
from reportlab.lib import colors
from reportlab.lib.pagesizes import letter, landscape
from reportlab.platypus import SimpleDocTemplate, Image, Spacer, Paragraph, PageBreak
from reportlab.platypus.tables import Table
from reportlab.lib.styles import getSampleStyleSheet
import datetime

path = "C:\\Users\\md1jdsp\\Desktop\\read_csv"

for filename in os.listdir(path):
    if filename.endswith(".csv"):
        file_path = os.path.join(path,filename)
        df = pd.read_csv(file_path)

output_filename = "C:\\Users\\md1jdsp\\Desktop\\read_csv\\output.pdf"
# Create a PDF document
doc = SimpleDocTemplate(output_filename, pagesize=landscape(letter))

# Create a list of flowables (elements) to add to the PDF
elements = []

# Create a title for the PDF report
# Title page content
title_style = getSampleStyleSheet()['Title']
title_text = "iBEAt study: results on the MR biomarkers"
subtitle_text = "Internal report"
subtitle_text_2 = datetime.datetime.now().strftime('%Y'+'/'+ '%m' +'/' +'%d')

# Title
elements.append(Paragraph(title_text, title_style))
elements.append(Spacer(1, 12))

# Subtitle
elements.append(Paragraph(subtitle_text, title_style))
elements.append(Spacer(1, 36))  # Add space after the subtitle

# Subtitle (date)
elements.append(Paragraph(subtitle_text_2, title_style))
elements.append(Spacer(1, 36))  # Add space after the subtitle

# Add a page breaker
elements.append(PageBreak())

# T1 # Generate and add a plot to the PDF
T1 = df[df["SeriesDescription"]=="T1"]
T1_BL = T1[T1["Baseline_Followup"]=="BL"]
T1_mean = T1_BL[T1_BL["Parameter"]=="Mean"]
T1_mean_LK = T1_mean[T1_mean["Region of Interest"]=="LK"]
T1_mean_RK = T1_mean[T1_mean["Region of Interest"]=="RK"]

plt.figure(figsize=(16, 9))

plt.title("iBEAt Whole Kidney Mean T1",fontsize=16)
plt.plot(T1_mean_LK["PatientID"],np.array(T1_mean_LK["Value"]),'>',label="Left Kidney",markersize=8)
plt.plot(T1_mean_RK["PatientID"],np.array(T1_mean_RK["Value"]),'<',label="Right Kidney",markersize=8)
plt.xticks(rotation = 90)
plt.yticks(fontsize = 10)
plt.ylim(np.min([np.array(T1_mean_LK["Value"]), np.array(T1_mean_RK["Value"])]), np.max([np.array(T1_mean_LK["Value"]), np.array(T1_mean_RK["Value"])]))
plt.rcParams['xtick.labelsize']=3

plt.axhline(y=1100, color='g', linestyle='-',linewidth=3)
plt.axhline(y=1400, color='g', linestyle='-',linewidth=3)
plt.legend()

plt.xlabel("Patient ID",fontsize=14)
plt.xticks(fontsize=4)
plt.ylabel("T1 [ms]",fontsize=14)

# Save the plot to a BytesIO buffer
plot_buffer = BytesIO()
plt.savefig(plot_buffer, format='png',dpi=600)
plt.close()

# Add the plot to the PDF as an Image
plot_image = Image(plot_buffer, width=800, height=450)
elements.append(plot_image)

# Add a page breaker
elements.append(PageBreak())

# T1 # Generate and add a plot to the PDF
T2 = df[df["SeriesDescription"]=="T2"]
T2_BL = T2[T2["Baseline_Followup"]=="BL"]
T2_mean = T2_BL[T2_BL["Parameter"]=="Mean"]
T2_mean_LK = T2_mean[T2_mean["Region of Interest"]=="LK"]
T2_mean_RK = T2_mean[T2_mean["Region of Interest"]=="RK"]

plt.figure(figsize=(16, 9))

plt.title("iBEAt Whole Kidney Mean T2",fontsize=16)
plt.plot(T2_mean_LK["PatientID"],np.array(T2_mean_LK["Value"]),'>',label="Left Kidney",markersize=8)
plt.plot(T2_mean_RK["PatientID"],np.array(T2_mean_RK["Value"]),'<',label="Right Kidney",markersize=8)
plt.xticks(rotation = 90)
plt.yticks(fontsize = 10)
plt.ylim(np.min([np.array(T2_mean_LK["Value"]), np.array(T2_mean_RK["Value"])]), np.max([np.array(T2_mean_LK["Value"]), np.array(T2_mean_RK["Value"])]))
plt.rcParams['xtick.labelsize']=3

plt.axhline(y=55, color='g', linestyle='-',linewidth=3)
plt.axhline(y=80 , color='g', linestyle='-',linewidth=3)
plt.legend()

plt.xlabel("Patient ID",fontsize=14)
plt.xticks(fontsize=4)
plt.ylabel("T2 [ms]",fontsize=14)

# Save the plot to a BytesIO buffer
plot_buffer = BytesIO()
plt.savefig(plot_buffer, format='png',dpi=600)
plt.close()

# Add the plot to the PDF as an Image
plot_image = Image(plot_buffer, width=800, height=450)
elements.append(plot_image)

# Add a page breaker
elements.append(PageBreak())

# T2s # Generate and add a plot to the PDF
T2s = df[df["SeriesDescription"]=="T2s"]
T2s_BL = T2s[T2s["Baseline_Followup"]=="BL"]
T2s_mean = T2s_BL[T2s_BL["Parameter"]=="Mean"]
T2s_mean_LK = T2s_mean[T2s_mean["Region of Interest"]=="LK"]
T2s_mean_RK = T2s_mean[T2s_mean["Region of Interest"]=="RK"]

plt.figure(figsize=(16, 9))

plt.title("iBEAt Whole Kidney Mean T2s",fontsize=16)
plt.plot(T2s_mean_LK["PatientID"],np.array(T2s_mean_LK["Value"]),'>',label="Left Kidney",markersize=8)
plt.plot(T2s_mean_RK["PatientID"],np.array(T2s_mean_RK["Value"]),'<',label="Right Kidney",markersize=8)
plt.xticks(rotation = 90)
plt.yticks(fontsize = 10)
plt.ylim(np.min([np.array(T2s_mean_LK["Value"]), np.array(T2s_mean_RK["Value"])]), np.max([np.array(T2s_mean_LK["Value"]), np.array(T2s_mean_RK["Value"])]))
plt.rcParams['xtick.labelsize']=3

plt.axhline(y=40, color='g', linestyle='-',linewidth=3)
plt.axhline(y=55 , color='g', linestyle='-',linewidth=3)
plt.legend()

plt.xlabel("Patient ID",fontsize=14)
plt.xticks(fontsize=4)
plt.ylabel("T2s [ms]",fontsize=14)

# Save the plot to a BytesIO buffer
plot_buffer = BytesIO()
plt.savefig(plot_buffer, format='png',dpi=600)
plt.close()

# Add the plot to the PDF as an Image
plot_image = Image(plot_buffer, width=800, height=450)
elements.append(plot_image)

# Add a page breaker
elements.append(PageBreak())

# FA # Generate and add a plot to the PDF
FA = df[df["SeriesDescription"]=="FA"]
FA_BL = FA[FA["Baseline_Followup"]=="BL"]
FA_mean = FA_BL[FA_BL["Parameter"]=="Mean"]
FA_mean_LK = FA_mean[FA_mean["Region of Interest"]=="LK"]
FA_mean_RK = FA_mean[FA_mean["Region of Interest"]=="RK"]

plt.figure(figsize=(16, 9))

plt.title("iBEAt Whole Kidney Mean FA",fontsize=16)
plt.plot(FA_mean_LK["PatientID"],np.array(FA_mean_LK["Value"]),'>',label="Left Kidney",markersize=8)
plt.plot(FA_mean_RK["PatientID"],np.array(FA_mean_RK["Value"]),'<',label="Right Kidney",markersize=8)
plt.xticks(rotation = 90)
plt.yticks(fontsize = 10)
plt.ylim(0.10, 0.35)
plt.rcParams['xtick.labelsize']=3

plt.axhline(y=0.15, color='g', linestyle='-',linewidth=3)
plt.axhline(y=0.30 , color='g', linestyle='-',linewidth=3)
plt.legend()

plt.xlabel("Patient ID",fontsize=14)
plt.xticks(fontsize=4)
plt.ylabel("FA [ms]",fontsize=14)

# Save the plot to a BytesIO buffer
plot_buffer = BytesIO()
plt.savefig(plot_buffer, format='png',dpi=600)
plt.close()

# Add the plot to the PDF as an Image
plot_image = Image(plot_buffer, width=800, height=450)
elements.append(plot_image)

# Add a page breaker
elements.append(PageBreak())

# ADC # Generate and add a plot to the PDF
ADC = df[df["SeriesDescription"]=="ADC"]
ADC_BL = ADC[ADC["Baseline_Followup"]=="BL"]
ADC_mean = ADC_BL[ADC_BL["Parameter"]=="Mean"]
ADC_mean_LK = ADC_mean[ADC_mean["Region of Interest"]=="LK"]
ADC_mean_RK = ADC_mean[ADC_mean["Region of Interest"]=="RK"]

plt.figure(figsize=(16, 9))

plt.title("iBEAt Whole Kidney Mean ADC",fontsize=16)
plt.plot(ADC_mean_LK["PatientID"],np.array(ADC_mean_LK["Value"]),'>',label="Left Kidney",markersize=8)
plt.plot(ADC_mean_RK["PatientID"],np.array(ADC_mean_RK["Value"]),'<',label="Right Kidney",markersize=8)
plt.xticks(rotation = 90)
plt.yticks(fontsize = 10)
plt.ylim(0.0010, 0.0023)
plt.rcParams['xtick.labelsize']=3

plt.axhline(y=0.0021, color='g', linestyle='-',linewidth=3)
plt.axhline(y=0.0015 , color='g', linestyle='-',linewidth=3)
plt.legend()

plt.xlabel("Patient ID",fontsize=14)
plt.xticks(fontsize=4)
plt.ylabel("ADC [ms]",fontsize=14)

# Save the plot to a BytesIO buffer
plot_buffer = BytesIO()
plt.savefig(plot_buffer, format='png',dpi=600)
plt.close()

# Add the plot to the PDF as an Image
plot_image = Image(plot_buffer, width=800, height=450)
elements.append(plot_image)

# Add a page breaker
elements.append(PageBreak())

# MT # Generate and add a plot to the PDF
MT = df[df["SeriesDescription"]=="MT"]
MT_BL = MT[MT["Baseline_Followup"]=="BL"]
MT_mean = MT_BL[MT_BL["Parameter"]=="Mean"]
MT_mean_LK = MT_mean[MT_mean["Region of Interest"]=="LK"]
MT_mean_RK = MT_mean[MT_mean["Region of Interest"]=="RK"]

plt.figure(figsize=(16, 9))

plt.title("iBEAt Whole Kidney Mean MT",fontsize=16)
plt.plot(MT_mean_LK["PatientID"],np.array(MT_mean_LK["Value"]),'>',label="Left Kidney",markersize=8)
plt.plot(MT_mean_RK["PatientID"],np.array(MT_mean_RK["Value"]),'<',label="Right Kidney",markersize=8)
plt.xticks(rotation = 90)
plt.yticks(fontsize = 10)
plt.ylim(np.min([np.array(MT_mean_LK["Value"]), np.array(MT_mean_RK["Value"])]), np.max([np.array(MT_mean_LK["Value"]), np.array(MT_mean_RK["Value"])]))
plt.rcParams['xtick.labelsize']=3

plt.axhline(y=35, color='g', linestyle='-',linewidth=3)
plt.axhline(y=25 , color='g', linestyle='-',linewidth=3)
plt.legend()

plt.xlabel("Patient ID",fontsize=14)
plt.xticks(fontsize=4)
plt.ylabel("MT [ms]",fontsize=14)

# Save the plot to a BytesIO buffer
plot_buffer = BytesIO()
plt.savefig(plot_buffer, format='png',dpi=600)
plt.close()

# Add the plot to the PDF as an Image
plot_image = Image(plot_buffer, width=800, height=450)
elements.append(plot_image)

# Add a page breaker
elements.append(PageBreak())

# RBF # Generate and add a plot to the PDF
RBF = df[df["SeriesDescription"]=="RBF"]
RBF_BL = RBF[RBF["Baseline_Followup"]=="BL"]
RBF_mean = RBF_BL[RBF_BL["Parameter"]=="Mean"]
RBF_mean_LK = RBF_mean[RBF_mean["Region of Interest"]=="LK"]
RBF_mean_RK = RBF_mean[RBF_mean["Region of Interest"]=="RK"]

plt.figure(figsize=(16, 9))

plt.title("iBEAt Whole Kidney Mean RBF (ASL)",fontsize=16)
plt.plot(RBF_mean_LK["PatientID"],np.array(RBF_mean_LK["Value"]),'>',label="Left Kidney",markersize=8)
plt.plot(RBF_mean_RK["PatientID"],np.array(RBF_mean_RK["Value"]),'<',label="Right Kidney",markersize=8)
plt.xticks(rotation = 90)
plt.yticks(fontsize = 10)
plt.ylim(np.min([np.array(RBF_mean_LK["Value"]), np.array(RBF_mean_RK["Value"])]), 350)
plt.rcParams['xtick.labelsize']=3

plt.axhline(y=100, color='g', linestyle='-',linewidth=3)
plt.axhline(y=300 , color='g', linestyle='-',linewidth=3)
plt.legend()

plt.xlabel("Patient ID",fontsize=14)
plt.xticks(fontsize=4)
plt.ylabel("RBF [ms]",fontsize=14)

# Save the plot to a BytesIO buffer
plot_buffer = BytesIO()
plt.savefig(plot_buffer, format='png',dpi=600)
plt.close()

# Add the plot to the PDF as an Image
plot_image = Image(plot_buffer, width=800, height=450)
elements.append(plot_image)

elements.append(PageBreak())

# FractionWater # Generate and add a plot to the PDF
FractionWater = df[df["SeriesDescription"]=="FrationWater"]
FractionWater_BL = FractionWater[FractionWater["Baseline_Followup"]=="BL"]
FractionWater_mean = FractionWater_BL[FractionWater_BL["Parameter"]=="Mean"]
FractionWater_mean_LK = FractionWater_mean[FractionWater_mean["Region of Interest"]=="LK"]
FractionWater_mean_RK = FractionWater_mean[FractionWater_mean["Region of Interest"]=="RK"]

plt.figure(figsize=(16, 9))

plt.title("iBEAt Whole Kidney Mean FractionWater",fontsize=16)
plt.plot(FractionWater_mean_LK["PatientID"],np.array(FractionWater_mean_LK["Value"]),'>',label="Left Kidney",markersize=8)
plt.plot(FractionWater_mean_RK["PatientID"],np.array(FractionWater_mean_RK["Value"]),'<',label="Right Kidney",markersize=8)
plt.xticks(rotation = 90)
plt.yticks(fontsize = 10)
plt.ylim(0.83,1)
plt.rcParams['xtick.labelsize']=3

plt.axhline(y=0.98, color='g', linestyle='-',linewidth=3)
plt.axhline(y=0.90 , color='g', linestyle='-',linewidth=3)
plt.legend()

plt.xlabel("Patient ID",fontsize=14)
plt.xticks(fontsize=4)
plt.ylabel("FractionWater [ms]",fontsize=14)

# Save the plot to a BytesIO buffer
plot_buffer = BytesIO()
plt.savefig(plot_buffer, format='png',dpi=600)
plt.close()

# Add the plot to the PDF as an Image
plot_image = Image(plot_buffer, width=800, height=450)
elements.append(plot_image)

doc.build(elements)
