#!/usr/bin/python3
import sys
import os
import math
import urllib.request
from requirements import *
from ftplib import FTP
from WaitingSpinnerWidget import QtWaitingSpinner
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *

def exec(cmd):
    return os.popen(cmd).read()

class Installer:
    def checkAptPackage(self, package):
        return exec("dpkg-query -W -f='${Status}' %s 2>/dev/null | grep -c 'ok installed'" % (package)).find('1') != -1
    
    def checkPipPackage(self, package, virtualenv=None):
        verSep = ['=', '>', '<']
        if package[:4] == 'git+':
            package = package[package.rfind('/')+1:]
        minSep = -1
        for i, c in enumerate(package):
            if c in verSep:
                minSep = i
        if minSep != -1:
            package = package[:minSep-1]

        if virtualenv == None:
            return len(exec("pip show %s" % (package))) != 0
        return len(exec('bash -c "source %s/bin/activate && pip show %s"' % (virtualenv, package))) != 0
        
    def installAptPackage(self, package):
        exec("sudo apt-get install -y %s" % (package))
        return self.checkAptPackage(package)

    def installPipPackage(self, package, virtualenv=None):
        if virtualenv == None:
            exec("pip install %s" % (package))
        else: exec('bash -c "source %s/bin/activate && pip install %s"' % (virtualenv, package))
        return self.checkPipPackage(package, virtualenv)



class InstallationThread(QThread):
    logSignal = pyqtSignal(str, bool)
    stopSignal = pyqtSignal()
    installer = Installer()
    def addLog(self, msg, overwriteLastWord=False):
        self.logSignal.emit(msg, overwriteLastWord)

    def downloadCB(self, count, blockSize, totalSize):
        act = count*blockSize
        if self.fsize != 0:
            self.addLog('%d%%(%.1fMB/%.1fMB)' % (int(act*100/self.fsize), act*0.000001, self.fsize*0.000001), True)
        else:
            if totalSize > 0:
                self.addLog('%d%%(%.1fMB/%.1fMB)' % (int(act*100/totalSize), act*0.000001, totalSize*0.000001), True)
            else:
                self.tsize += blockSize
                self.addLog('%dBytes' % int(self.tsize), True) 

######################################################
#################INSTALLATION ROUTINE#################
######################################################
    def run(self):
        for p in systemPkgs:
            self.addLog('Installing system package: %s...' % (p))
            if not self.installer.installAptPackage(p):
                self.stopSignal.emit()
                self.addLog('Error installing system package: %s' % (p))
                return
            

        self.addLog('Creating Python virtual environment...')
        exec('virtualenv --python=/usr/bin/%s ../env' % (pythonVer))


        for p in pipPkgs:
            self.addLog('Installing pip package: %s...' % (p))
            if not self.installer.installPipPackage(p, '../env'):
                self.stopSignal.emit()
                self.addLog('Error installing pip package: %s' % (p))
                return
            

        self.addLog('Cloning peptide repo...')
        if os.path.exists('../peptide'):
            exec('rm -rf ../peptide')
        exec('git clone https://github.com/mchelem/peptide ../peptide')
        self.addLog('Compiling torsions lib...')
        exec('cc -o ../libs/torsions ../libs/torsions.c -lm')
        if not os.path.isfile('../libs/torsions'):
            self.stopSignal.emit()
            self.addLog('Error compiling torsions library!')
            return


        self.addLog('Downloading PDB SEQRES database: 0%')
        exec('mkdir -p ../data/blastdb')#ftp://ftp.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt.gz
        ftp = FTP('ftp.wwpdb.org')
        ftp.login()
        self.fsize = ftp.size('/pub/pdb/derived_data/pdb_seqres.txt.gz')
        ftp.close()
        urllib.request.urlretrieve('ftp://ftp.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt.gz', \
            '../data/blastdb/pdb_seqres.txt.gz', reporthook=self.downloadCB)
        self.addLog('Building BLAST database...')
        exec('cd ../data/blastdb && gunzip -f pdb_seqres.txt.gz && makeblastdb -in pdb_seqres.txt -out pdbseqres -dbtype prot')
        self.fsize = 0


        self.addLog('Downloading precomputed databases (Secondary Structure): 0%')
        urllib.request.urlretrieve('https://www.dropbox.com/s/xb7qc49p75k9qlf/ss.db.tgz?dl=1', \
            '../ss.db.tgz', reporthook=self.downloadCB)
        self.addLog('Downloading precomputed databases (Torsions): 0%')
        urllib.request.urlretrieve('https://www.dropbox.com/s/18q2t35eh6rml1s/torsions.db.tgz?dl=1', \
            '../torsions.db.tgz', reporthook=self.downloadCB)     
        self.addLog('Extracting precomputed structures databases...')
        exec('tar xzvf ../ss.db.tgz -C ../')
        exec('rm -f ../ss.db.tgz')
        exec('tar xzvf ../torsions.db.tgz -C ../')
        exec('rm -f ../torsions.db.tgz')


        self.addLog('Downloading SCRATCH-1D: 0%')
        urllib.request.urlretrieve('http://download.igb.uci.edu/SCRATCH-1D_1.0.tar.gz', \
            '../libs/SCRATCH-1D_1.0.tar.gz', reporthook=self.downloadCB)
        self.addLog('Installing SCRATCH-1D... (may take a while)')
        exec('cd ../libs && tar xvzf SCRATCH-1D_1.0.tar.gz && cd SCRATCH-1D_1.0 && perl install.pl')
        self.addLog('---------------------------------------------------')
        self.addLog('CReF Installed!')
        self.addLog("Run './cref.sh --help' for more information.")
        self.stopSignal.emit()
        pass
######################################################
######################################################
######################################################  

class GUI(QWidget):
    def __init__(self):
        super().__init__()
        exec('sudo echo')
        self.initGUI()
        self.fsize = 0
        self.tsize = 0

    def addLog(self, msg, overwrite):
        if overwrite:
            t = self.logEditor.toPlainText()
            t = t[:t.rfind(' ')] + ' ' + msg
            self.logEditor.setPlainText(t)
        else:
            self.logEditor.append(msg)
        self.logEditor.verticalScrollBar().setValue(
            self.logEditor.verticalScrollBar().maximum()
        )

    def stop(self):
        self.loading.stop()
        self.loading.setVisible(False)


    def installBtnCB(self):
        self.install = InstallationThread()
        self.install.start()
        self.install.logSignal.connect(self.addLog)
        self.install.stopSignal.connect(self.stop)
        self.installBtn.setVisible(False)
        self.loading.setVisible(True)
        self.loading.start()
        

   
    def initGUI(self):           
        self.resize(483, 579)
        fg = self.frameGeometry()
        center = QDesktopWidget().availableGeometry().center()
        fg.moveCenter(center)
        self.move(fg.topLeft())
        self.setWindowTitle('CReF2')   
        

        self.installBtn = QPushButton('Install', self)
        self.installBtn.move(202, 545)
        self.installBtn.clicked.connect(self.installBtnCB)

        self.logEditor = QTextEdit(self)
        self.logEditor.setReadOnly(True)
        self.logEditor.resize(463, 206)
        self.logEditor.move(10, 330)

        logo = QLabel(self)
        pm = QPixmap('logo.png')
        logo.setPixmap(pm)
        logo.move(10, 10)

        self.loading = QtWaitingSpinner(self, False, False)
        self.loading.setRoundness(70.0)
        self.loading.setMinimumTrailOpacity(15.0)
        self.loading.setTrailFadePercentage(70.0)
        self.loading.setNumberOfLines(30)
        self.loading.setLineLength(2)
        self.loading.setLineWidth(5)
        self.loading.setInnerRadius(10)
        self.loading.setRevolutionsPerSecond(1)
        self.loading.setColor(QColor(42, 161, 95))
        self.loading.move(self.width() / 2 - self.loading.width() / 2,545)
        self.loading.setVisible(False)
                

        self.show()
        
if __name__ == '__main__':
    app = QApplication(sys.argv)
    w = GUI()
    sys.exit(app.exec_())