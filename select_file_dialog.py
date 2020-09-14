from sys import executable, argv
from subprocess import check_output
from PyQt5.QtWidgets import QFileDialog, QApplication

def gui_fname(directory='./', filter="iqa"):
    """Open a file dialog, starting in the given directory, and return
    the chosen filename"""
    if filter == "iqa":
        filter = "JSON (*.json);;GZIPPED JSON (*.json.gz);;MessagePack (*.iqapack);;All files (*)"
    elif filter == "surf":
      filter = "Surface file (*.h5);;All files (*)"
    elif filter == "dgrid_wfn":
      filter = "FHI-aims special files (*.fhi);;ADF TAPE21 (*.adf);;ADF TAPE21 loc. orbitals (*.adf_loc);;GAMESS output (*.gms);;GAUSSIAN94 fchkpt (*.g94);;GAUSSIAN98 fchkpt (*.g98);;GAUSSIAN03 fchkpt (*.g03);;GAUSSIAN fchkpt (*.g09);;GAUSSIAN WFN (*.gwf);;MOLPRO molden file (*.mp);;MOLCAS molden file (*.mc);;MOLDEN file (*.md);;CLEMENTI-ROETTI tables (*.CR);;All files (*)"

    # run this exact file in a separate process, and grab the result
    file = check_output([executable, __file__, directory, filter])
    return file.strip()

if __name__ == "__main__":
    directory = argv[1]
    filter   = argv[2]
    app = QApplication([directory])
    fname = QFileDialog.getOpenFileName(None, "Select a file...", 
            directory, filter=filter)
    print(fname[0])
