# PyLink
### a PyMOL plugin to identify links

VERSION: 1.3 
Compatible with Python 3 and PyMOL 2.x
Updated on 25.05.2021

The PYMOL PyLink plugin requires the following Python packages. 
Make sure they are installed before installing this plugin in pymol:

    numpy
    matplotlib
    pmw

These packages may be installed using the following commands.

    pip3 install -r requirements.txt

Details @ [pylink.cent.uw.edu.pl](http://pylink.cent.uw.edu.pl)

### Install plugin in PyMOL

Open PyMOL and go to: 

plugin -> plugin manager -> install new plugin -> choose file 

Select the 

    __init__.py 

file from the decompressed folder containing the downloaded PyLasso plugin.


### Problem with outdated Python pmw package

If you see the following error after opening PyMOL:

    Traceback (innermost last):
      File "C:\Python38\lib\site-packages\Pmw\Pmw_2_0_1\lib\PmwBase.py", line 1776, in __call__
        return self.func(*args)
      File "C:\Python38\lib\site-packages\pymol\plugins\legacysupport.py", line 85, in plugin_manager
        managergui.manager_dialog()
      File "C:\Python38\lib\site-packages\pymol\plugins\managergui.py", line 43, in manager_dialog
        dialog = PluginManager(get_tk_root())
      File "C:\Python38\lib\site-packages\pymol\plugins\managergui.py", line 125, in __init__
        notebook = Pmw.NoteBook(master)
      File "C:\Python38\lib\site-packages\Pmw\Pmw_2_0_1\lib\PmwNoteBook.py", line 60, in __init__
        Pmw.Color.bordercolors(self, self['hull_background'])
      File "C:\Python38\lib\site-packages\Pmw\Pmw_2_0_1\lib\PmwColor.py", line 359, in bordercolors
        '#%04x%04x%04x' % (lightRGB[0], lightRGB[1], lightRGB[2]),
    TypeError: %x format: an integer is required, not float

that means you have an old version of the Python PMW package. You need to install the patched one from:

https://github.com/schrodinger/pmw-patched


## PyMOL on Mac OS:

This setup was tested on Mac OS X 11.4 (Big Sur)

Install Xquartz from: https://www.xquartz.org/

Install PyMOL using MacPorts (https://www.macports.org/)

    sudo port install pymol py39-matplotlib py39-pmw

After installation close the terminal and open a new one
Then start PyMOL by typing:
    
    pyymol

