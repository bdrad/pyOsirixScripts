# A DICOM-Embedded Annotation System for 3D Cross-Sectional Imaging Data
A Heydari, BA, Berkeley, CA; D Lituiev; T H Vu, MD, PhD; Y Seo, PhD; J Sohn, MD (aarashy@gmail.com)
 
## BACKGROUND
Annotation of imaging data, especially in 3D cross sectional imaging, remains challenging due to a paucity of common storage standards, user-friendly annotation tools, and visualization methods. Programmatic visualization tools such as Matplotlib and unspecialized formats such as CSV and XML are unfavorable for visualizing 3D cross-sectional annotations. We present a data-generation workflow that embeds radiologist annotations directly into DICOM headers. The resulting annotations are ready for use in machine learning projects and widely compatible and easily visualized in DICOM viewers such as OsiriX/Horos.

 
## EVALUATION
Our python-based annotation pipeline leverages the DICOM header at tag value (60xx, 3000), which supports 16 binary mask ROIs per cross-sectional image. It extracts radiologist annotations from OsiriX/Horos via the PyOsiriX toolkit and embeds them directly into the DICOM header. Additionally, a reverse pipeline was developed to push annotations from CSV databases into DICOM headers to enable visualization within DICOM viewers. This process is demonstrated with lung nodule segmentations in chest CTs in Horos (Figure).
 
## DISCUSSION
Effective annotation and visualization are among the most important factors towards successful deep learning model development. Our tool supports all annotation shapes, including oval, rectangular, polygon, and penciled ROIs. Leveraging the DICOM header overlay tags makes data visualization intuitive and compatible with most DICOM viewers. Our pipelines bridge between radiologist annotations in DICOM viewers and machine learning friendly formats, ultimately facilitating radiological machine learning research.

## CONCLUSION
A widely compatible DICOM-native annotation pipeline with strong support for data visualization was created to facilitate computer vision projects involving 3D cross-sectional data.
