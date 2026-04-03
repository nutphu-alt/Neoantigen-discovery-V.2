# Neoantigen-discovery-V.2

## This project aims to create an automatic pipeline for neoantigen discovery. It is under development.  

### This pipline is adapted from Kreiter, Vormehr et al. 2015.

The overall pipeline is illustrated below. It begins with performing RNA and exome sequencing from tumor tissue or cancer cell line, in order to identify mutations. Then, exome sequencing is performed on normal tissue such as human PBMC or mouse tail tip. Mutations found in normal tissue are eliminated from mutations found on tumor tissue or cancer celline to substract germline mutations.
<br>

<p align="center">
<img width="600" height="300" alt="image" src="https://github.com/user-attachments/assets/8e8847ec-85fc-4633-ba57-d394bb87c10c" />  
</p>
  
The steps, used to indentify neoantigen, are listed below.
<br>
<br>

<p align="center">
<img width="300" height="750" alt="image" src="https://github.com/user-attachments/assets/4f722fd9-f7e6-497f-9369-c0f501d72936" />
</p>
