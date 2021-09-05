# stat-thermo

## How to run:

### 1. Download the files from github
    
    Click on the green "Code" button on this page.
    
    Select "Download .zip"
    
    Unarchive the downloaded .zip

### 2. Make sure you have python3
  In Powershell (Windows) or Bash (Linux) type: python --version
  
  If it isn't installed then install it:
  
    Windows: https://www.microsoft.com/fr-fr/p/python-37/9nj46sx7x90p?rtc=1&activetab=pivot:overviewtab
    
    Linux (Ubuntu): Type in bash: sudo apt-get update
                                  sudo apt-get install python3
                                  
    Mac should come with python out-of-the-box. 

### 3. Open Powershell (or terminal) and run the command below:
    
    To open Powershell:
        Windows + R
        Type "powershell" and hit enter
        cd + path of the directory you have main.py saved in
            example: cd Desktop

    python3 main.py HF.txt HF2.txt
  
### 4. For the above command, you can replace HF.txt and HF2.txt with other files

        example: python3 main.py H2.txt CO.txt H2CO.txt

### 5. The data in the .txt is as follows:

            Molecule name
            Reaction coefficient (has to be an integer!)
            Liniarity (type linear if linear, non-linear otherwise)
            Sigma
            Mass of molecule
            Rotational contsants (B) values in a row, split by space, in cm^-1 (if N/A, type 0)
            Vibration wavenumbers values in a row, split by space, in cm^-1 (if N/A, type 0)
            Vibrational degeneracies in a row, split by space (if N/A, type 0)
            (Optional) Electronic wavenumbers values in a row, split by space, in cm^-1
            (Optional) Electronic degeneracies in a row, split by space

### 6. Follow instructions on screen and enjoy! :D
