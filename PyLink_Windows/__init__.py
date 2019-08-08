# -*- coding: utf-8 -*-
# PyLink: A PyMOL plugin to identify link, version for Windows platform
# Script/plugin by Aleksandra Gierut (a.gierut@cent.uw.edu.pl)
# Questions should be addressed to: Joanna Sulkowska (jsulkowska@cent.uw.edu.pl)
# Any technical difficulties and remarks please report to: Aleksandra Gierut (a.gierut@cent.uw.edu.pl)
#
# THE AUTHOR(S) DISCLAIM ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
# INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS.  IN
# NO EVENT SHALL THE AUTHOR(S) BE LIABLE FOR ANY SPECIAL, INDIRECT OR
# CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF
# USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR
# OTHER TORTUOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
# PERFORMANCE OF THIS SOFTWARE.
# ----------------------------------------------------------------------
import re
import shutil
import textwrap
import random

from itertools import combinations
from math import floor
import glob
from tkFont import Font
import os

plugin_path = os.path.dirname(__file__)
system_working_directory = os.getcwd()

try:
    import Pmw
    import Tkinter as tk
    import tkFileDialog
except:
    print "  ### Graphic libraries not found. Please install them (Tkinter and Pmw) and re-run the plugin."

from pymol.cgo import *
from pymol import cmd
from pymol import util

try:
    import matplotlib as mplt
    mplt.use('TkAgg')
    import matplotlib.pyplot as pyplt
    import matplotlib.image as mpimg
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
    from matplotlib.lines import Line2D
except:
    print "  ### Matplotlib library not found. Please install it and re-run the plugin."


def __init__(self):
    self.menuBar.addmenuitem('Plugin', 'command', 'PyLink',
                             label='PyLink',
                             command=lambda s=self: PyLink(s))


class PyLink:

    def __init__(self, app):
        self.parent = app.root
        if cmd.get_version()[1] < 2.0:
            self.bold_font = Font(family="Helvetica", size=8, weight="bold")
            self.row_name_font = "bold 10"
            self.surf_btn_width = 15
            self.gln_fig_dpi = 98
            self.gln_fig_size = (7.8, 3.5)
            self.chart_img_axes = [0.3, 0.4, 0.4, 0.35]
            self.chart_size = 3.3
            self.chart_font_size = 28
            self.chart_subtitle = 16
        else:
            self.parent.option_add("*Font", "Helvetica 9")
            self.bold_font = Font(family="Helvetica", size=8, weight="bold")
            self.row_name_font = "Helvetica 9"
            self.surf_btn_width = 18
            self.gln_fig_dpi = 98
            self.gln_fig_size = (9, 4)
            self.chart_img_axes = [0.3, 0.4, 0.42, 0.38]
            self.chart_size = 3.0
            self.chart_font_size = 18
            self.chart_subtitle = 14

        # Link variables
        self.screen_height = self.parent.winfo_screenheight()
        self.screen_width = self.parent.winfo_screenwidth()
        self._filenames = []
        self._chains = []
        self._is_file_xyz = {}
        self._is_file_macrolink = {}
        self._contains_bridge = {}
        self._chains_marginal_atoms = {}
        self._marginal_atoms = {}
        self.is_macrolink = False
        self.is_ions = False
        self.prog_converter = "python " + plugin_path + os.sep + 'converter.py '
        self.prog_homfly = plugin_path + os.sep + "homfly.exe "
        self.prog_poly = plugin_path + os.sep + "poly.exe  "
        self.prog_ncuc = plugin_path + os.sep + "ncucLinks.exe "
        self.prog_gln = "python " + plugin_path + os.sep + "GLNtoPNG.py "
        self.view_btn_width = 3
        self.distance_error = {}
        self.list_hint_distance = []
        self.previously_displayed_bond = []
        self.bottom_btns_width = 6
        self.suggest_btn_width = 18
        self.retrieve_btns_width = 9
        self.retrieve_btns_pad = 1 # needed?
        self.hint_width = 60
        # Link information window
        self.protein_row_width = 40
        self.protein_row_name_width = 10
        self.protein_btn_width = 20
        self.view_details_btn_width = 10
        self.gln_row_name_width = 28
        # Advanced variables
        self.is_correct_data = tk.IntVar()
        self.is_gln_selected = tk.IntVar()
        self.whole_chains_enabled = tk.IntVar()
        # Macromolecular links
        self.joining_method = tk.IntVar()

        self.xyz_chain_possibilities = list(reversed((string.ascii_uppercase))) # 2+ XYZ files needs to have chains assigned

        self.colors_of_piercings = [[0.0, 0.0, 1.0], [1.0, 1.0, 0.0], [1.0, 1.0, 0.0], [1.0, 0.0, 0.0]]
        self.colors_of_sequence = ["yelloworange", "violet", "lightteal", "limegreen", ]
        self.colors_of_surfaces = [[1.0, 1.0, 0.5], [1.0, 0.75, 0.87], [0.8, 1.0, 1.0], [0.65, 0.9, 0.65]]

        self.load_file()
        if self._filenames:
            self.initialise_plugin_interface()
            self.adjust_object_representation()

    ####################################################################################################################
    #                             CHECK FILE EXTENSION & ADJUST POLYMER REPRESENTATION
    ####################################################################################################################

    def load_file(self):
        self.open_file_window = tkFileDialog.askopenfilenames(initialdir=os.getcwd(), title="PyLink",
                                                              filetypes=[("PDB", "*.pdb"), ("XYZ", "*.xyz")])

        self._full_path_to_files = []
        if not self.open_file_window or None == self.open_file_window:
            print "# No file was chosen. PyLink was shut down."
            return
        elif len(self.open_file_window) > 4:
            print "# Too many files were chosen. Please pick at most four XYZ files and one PDB file."
            return
        elif all([path_file.endswith(".pdb") for path_file in self.open_file_window]) and \
                [path_file.endswith(".pdb") for path_file in self.open_file_window].count(True) > 1:
            print "# Too many PDB files were chosen. Choose only one file and try again."
            return
        elif 1 < len(self.open_file_window) <= 4:
            for i in self.open_file_window:
                i = i.replace("/", "\\")
                self._full_path_to_files.append(i)
        else:
            self._full_path_to_files.append(self.open_file_window[0].replace("/", "\\"))
        self._full_path_to_dir = os.sep.join(self._full_path_to_files[0].split(os.sep)[:-1])

        cmd.reinitialize()
        for idx, filepath in enumerate(self._full_path_to_files):
            filename = filepath.split(os.sep)[-1]
            if filename.endswith("pdb"):
                if self.is_input_marcolink(idx):
                    self.is_macrolink = True
                    filename = filepath.split(os.sep)[-1][:-4] + "_CA.pdb"
                    self.simplify_macrolink_to_ca(filepath)
                    self._filenames.append(filename)
                    self._full_path_to_files[idx] = self._full_path_to_dir + os.sep + self._filenames[-1]
                    self._is_file_xyz[filename[:-4]] = False
                    self._is_file_macrolink[filename[:-4]] = True
                else:
                    self._filenames.append(filepath.split(os.sep)[-1])
                    self._is_file_xyz[filename[:-4]] = False
                    self._is_file_macrolink[filename[:-4]] = False
                self._contains_bridge[filename[:-4]] = self.contains_bridge_data(idx)
            else:
                self._filenames.append(filepath.split(os.sep)[-1][:-4] + "_xyz2.pdb")
                self._is_file_xyz[self._filenames[-1][:-4]] = True
                self._is_file_macrolink[self._filenames[-1][:-4]] = False
                self._contains_bridge[self._filenames[-1][:-4]] = False
                self.convert_xyz_to_pdb(self._full_path_to_files[idx], self._full_path_to_dir, self._filenames[idx])
                self._full_path_to_files[idx] = self._full_path_to_dir + os.sep + self._filenames[-1]

            cmd.load(filename=self._full_path_to_files[idx])
            structure = []
            for i in cmd.get_chains(self._filenames[idx][:-4]):
                structure.append(i)
            self._chains.append(structure)
        print "  Structures reloaded..."

    def is_input_marcolink(self, file_idx):
        input_file = open(self._full_path_to_files[file_idx], "r").read()
        expdata = list(re.findall('^MODEL|SOLUTION NMR', input_file, flags=re.M | re.S))

        return "MODEL" in expdata and "SOLUTION NMR" not in expdata

    def simplify_macrolink_to_ca(self, filename):
        macrolink_ca_only = open(filename[:-4] + "_CA.pdb", "w")

        with open(filename, "r") as protein:
            data = re.compile("^(SSBOND|LINK|TER|ATOM.*\d  CA|MODEL|ENDMDL)")
            for line in protein:
                if data.match(line):
                    macrolink_ca_only.write(line)
        macrolink_ca_only.close()

    def contains_bridge_data(self, file_idx):
        input_file = open(self._full_path_to_files[file_idx], "r").read()
        ss_bonds = list(re.findall('SSBOND|LINK', input_file, flags=re.M | re.S))

        return len(ss_bonds) is not 0

    def convert_xyz_to_pdb(self, current_file, path_to_xyz_file, new_xyz_file, ch=None):
        output_file = open(path_to_xyz_file + os.sep + new_xyz_file, 'w')
        chain = self.xyz_chain_possibilities.pop() if not ch else ch

        idx = 0
        xyz_atm = re.compile(
            r"^\s*(?P<resid>[0-9]*)\s+(?P<x>-?\d+\.\d*)\s+(?P<y>-?\d+\.\d*)\s+(?P<z>-?\d+\.\d*).*$")
        with open(current_file) as f:
            for line in f:
                data = xyz_atm.match(line)
                if data:
                    idx += 1
                    resid = int(data.groups()[0])
                    x = float(data.groups()[1])
                    y = float(data.groups()[2])
                    z = float(data.groups()[3])
                    output_file.write("ATOM%7d  CA  GLY %s%4d    %8.3f%8.3f%8.3f  1.00  0.00           C\n" % (
                        resid, chain, resid, x, y, z))
            output_file.write("TER     %d      GLY %s %d" % (resid + 1, chain, resid))
        output_file.close()

    def adjust_object_representation(self):
        for idx, link_file in enumerate(self._full_path_to_files):
            self.retrieve_marginal_atoms(link_file)
            file = self._filenames[idx][:-4]

            if not self._is_file_macrolink[file] and (self._is_file_xyz[file] or not self._contains_bridge[file]):
                self.connect_xyz_points(file)
                cmd.show(representation="lines", selection=file)
                cmd.set(name="line_width", value="4")
            else:
                cmd.hide(representation="nonbonded", selection="solvent")
                cmd.hide(representation="lines", selection=file)
                cmd.show(representation="cartoon", selection=file)
                cmd.cartoon(type="tube", selection=file)
            if self._is_file_macrolink[file]:
                cmd.hide(representation="cartoon", selection="all")
                cmd.show(representation="nonbonded", selection="all")
                cmd.split_states(object=file)
                util.cbc(selection="all")
                cmd.center(selection=file)
                cmd.orient(selection=file)
                cmd.clip(mode="slab", distance="1000")

        cmd.set(name="stick_color", value="orange", selection="all")
        cmd.set(name="sphere_color", value="orange", selection="all")
        cmd.set(name="sphere_scale", value=0.7)
        cmd.set(name="stick_radius", value=0.5)
        util.cbc()

    def retrieve_marginal_atoms(self, path_to_file):
        atoms = []
        reg = re.compile('ATOM\s\s+\d+')
        chain = ""
        beg = ""
        link_file = path_to_file.split(os.sep)[-1][:-4]

        with open(path_to_file, "r") as f:
            for line in f:
                clear_line = filter(len, line.split(" "))
                if clear_line[0] == "TER":
                    if clear_line[4].isdigit():
                        self._chains_marginal_atoms[link_file + chain] = (beg, clear_line[4])
                    else:
                        self._chains_marginal_atoms[link_file + chain] = (beg, clear_line[3][1:])
                    chain = ""
                if reg.match(line):
                    if clear_line[4].isalpha():
                        atoms.append(clear_line[5])
                        if chain == "":
                            chain = clear_line[4]
                            beg = clear_line[5]
                    else:
                        idx = clear_line[4][1:]
                        atoms.append(idx)
                        if chain == "":
                            chain = clear_line[4][0]
                            beg = clear_line[4][1:]
        self._marginal_atoms[link_file] = ([int(atoms[0]), int(atoms[-1])])

    def connect_xyz_points(self, filename):
        for i in range(self._marginal_atoms[filename][0], self._marginal_atoms[filename][1]):
            cmd.bond(atom1=filename + " and id " + str(i) + " and name ca",
                     atom2=filename + " and id " + str(i + 1) + " and name ca")

    ####################################################################################################################
    #                                            INITIALISE INTERFACE
    ####################################################################################################################

    def initialise_plugin_interface(self):
        link_title = ''
        for i in self._filenames:
            link_title += i + ', '
        link_title = '[' + link_title[:-2] + ']'

        self.dialog = Pmw.Dialog(self.parent, buttons=('Proceed', 'Exit'), title='PyLink ' + link_title + '',
                                 command=self._invoke_plugin_action)
        self.dialog.geometry("+%d+%d" % (self.screen_width / 4, self.screen_height / 4))

        self.notebook = Pmw.NoteBook(self.dialog.interior(), raisecommand=lambda x: self.create_notebook_page(x))
        self.notebook.add("Deterministic")
        self.notebook.add("Probabilistic")
        self.notebook.add("Macromolecular links")

        if all(self._is_file_macrolink[struct] for struct in self._is_file_macrolink.keys()):
            self.notebook.selectpage("Macromolecular links")
            self.disable_pages = True
        else:
            self.disable_pages = False
            cmd.set(name="seq_view", value=1)
        cmd.set(name="mouse_selection_mode", value=6)
        self.notebook.pack(fill='both', expand=1, padx=5, pady=5)

        self.dialog.resizable(0, 0)
        self.dialog.show()

    def create_notebook_page(self, page):
        self.type_closing_choice = [False, False]

        if page == "Deterministic":
            self.create_deterministic_page()
        elif page == "Probabilistic":
            self.create_probabilistic_page()
        else:
            self.create_macromolecular_page()

    def create_deterministic_page(self):
        if hasattr(self, "disable_pages") and self.disable_pages:
            self.notebook.selectpage("Macromolecular links")
            self.raise_popup_menu("No deterministic approach for macromolecular links.")

        self.hide_interface_components()

        self.page_deterministic = tk.Label(self.notebook.page(0))
        self.page_deterministic.grid(sticky='eswn', column=0, row=0)

        self.group_method = Pmw.Group(self.page_deterministic, tag_text='Link identification method')
        self.group_method.grid(sticky='eswn', column=0, row=0, padx=5, pady=5)

        self.identification_method = Pmw.ComboBox(self.group_method.interior(), entry_width=30,
                                                  scrolledlist_items=('automatic detection',
                                                                      'choose bridging atoms',
                                                                      'choose bridging atoms (extended on ions)'),
                                                  entryfield_entry_state="readonly",
                                                  selectioncommand=lambda x: self.set_link_closing())
        self.identification_method.selectitem(0)
        self.identification_method.grid(sticky='eswn', column=0, row=1, padx=5, pady=5)

        self.group_advanced = Pmw.Group(self.page_deterministic, tag_text='Advanced options')
        self.group_advanced.grid(sticky='eswn', column=1, row=0, padx=5, pady=5, rowspan=3)

        self.create_smooth_interior()
        self.create_structure_correctness_interior()
        self.create_gln_interior()

        hint_advanced = Pmw.Balloon(self.page_deterministic, relmouse="both")
        self.create_advanced_frame_hints(hint_advanced)

        self.notebook.setnaturalsize()

    def raise_popup_menu(self, error_message):
        if hasattr(self, "error_pop_menu") and self.error_popup.winfo_exists():
            self.error_popup.withdraw()

        self.error_popup = Pmw.MessageDialog(self.parent, title='Error!', defaultbutton=0,
                                             message_text=textwrap.fill(error_message, 80))
        self.error_popup.geometry("+%d+%d" % (self.screen_width / 2 - 150, self.screen_height / 2))

        self.error_popup.focus_force()
        self.error_popup.wait_window()

    def hide_interface_components(self):
        if hasattr(self, "chosen_own_link"):
            self.chosen_own_link.grid_forget()
        if hasattr(self, "chosen_own_ions_link"):
            self.chosen_own_ions_link.grid_forget()
        if hasattr(self, "chosen_automatic_macrolink"):
            self.chosen_automatic_macrolink.grid_forget()
        if hasattr(self, "chosen_automatic_macrolink"):
            self.chosen_automatic_macrolink.grid_forget()
        if hasattr(self, "group_number_of_projections") and self.notebook.getcurselection != "Probabilistic":
            self.group_number_of_projections.grid_forget()

    def set_link_closing(self):
        if str(self.identification_method.get()) == "automatic detection":
            self.choose_automatic_closing()
        elif str(self.identification_method.get()) == "choose bridging atoms":
            self.choose_own_closing()
        else:
            self.choose_ions_closing()
        self.notebook.setnaturalsize()

    def choose_automatic_closing(self):
        if hasattr(self, "chosen_own_link"):
            self.chosen_own_link.grid_forget()
            self.type_closing_choice[0] = False
        if hasattr(self, "chosen_own_ions_link"):
            self.chosen_own_ions_link.grid_forget()
            self.type_closing_choice[1] = False
        self.is_ions = False

    def choose_own_closing(self):
        if hasattr(self, "chosen_own_ions_link"):
            self.chosen_own_ions_link.grid_forget()
            self.type_closing_choice[1] = False
        if self.type_closing_choice[0] == False:
            self.type_closing_choice[0] = True
            self.is_ions = False
            cmd.set(name="mouse_selection_mode", value=6)
            cmd.set(name="seq_view", value=1)

            self.chosen_own_link = tk.Frame(self.group_method.interior())
            self.chosen_own_link.grid(sticky='swen', column=0, row=2, padx=5, pady=5)

            loops_list_label = tk.Label(self.chosen_own_link, text=textwrap.fill('Define the indices of residues '
                                                                                 'closing the loop. Can pick the part '
                                                                                 'or whole chains in this checkbuttons '
                                                                                 'next to the end of the chain.', 70))
            loops_list_label.grid(column=0, row=0, columnspan=5, pady=2, padx=5)

            tmp_frame = tk.Frame(self.chosen_own_link, padx=0, pady=0)
            self.btns_get_data = Pmw.ButtonBox(tmp_frame, orient="vertical", pady=self.retrieve_btns_pad)
            self.btns_get_data.add("Get from\nstructure", width=self.retrieve_btns_width, font=self.bold_font,
                                   command=lambda: self.get_loop_from_protein_structure())
            self.btns_get_data.add("Get from\nsequence", width=self.retrieve_btns_width, font=self.bold_font,
                                   command=lambda: self.get_loop_from_protein_sequence())
            self.btns_get_data.grid(sticky='wen', column=0, row=0)
            tmp_frame.grid(column=0, row=2, rowspan=4)

            chain_tuple = []
            for idx, elem in enumerate(self._filenames):
                for i in self._chains[idx]:
                    chain_tuple.append(elem[:-4] + " " + i)

            self.loops_list = []
            self.btns_view = []
            for j in range(4):
                self.loops_list.append((
                    (Pmw.ComboBox(self.chosen_own_link, labelpos='w', label_text='Chain', entry_width=10,
                                  scrolledlist_items=tuple(chain_tuple),
                                  selectioncommand=lambda x=j: self.pymol_display_chain([x]))),
                    (Pmw.EntryField(self.chosen_own_link, entry_width=7, entry_state=tk.DISABLED,
                                    validate={'validator': 'integer'})),
                    (Pmw.EntryField(self.chosen_own_link, entry_width=7, entry_state=tk.DISABLED,
                                    validate={'validator': 'integer'}))
                ))
                self.loops_list[-1][0].grid(column=1, row=j + 3, pady=10)
                self.loops_list[-1][0].configure(entryfield_entry_state="readonly")
                self.loops_list[-1][0].component("entryfield").insert(0, chain_tuple[0])
                self.loops_list[-1][1].grid(column=2, row=j + 3, pady=10)
                self.loops_list[-1][2].grid(column=3, row=j + 3, pady=10)

                self.btns_view.append(tk.Button(self.chosen_own_link, text="View", width=self.view_btn_width,
                                                state=tk.DISABLED,
                                                command=lambda x=self.loops_list[-1], y=j: self.view_possible_loop(x, y)))
                self.btns_view[-1].grid(sticky='w', column=4, row=j + 3)

            self.btn_whole_chains = tk.Checkbutton(self.chosen_own_link, variable=self.whole_chains_enabled,
                                                   text="Use whole chains",
                                                   command=lambda: self.modify_whole_chains())
            self.btn_whole_chains.grid(sticky='w', column=2, row=2, columnspan=3)
            self.btn_whole_chains.select()

            self.btns_bottom = Pmw.ButtonBox(self.chosen_own_link, orient="horizontal")
            self.btns_bottom.add("Clear", width=self.bottom_btns_width, command=lambda: self.clear_loops_list())
            self.btns_bottom.grid(sticky='swn', column=2, row=7, columnspan=5, padx=5, pady=2)

            hint_choose_own_link = Pmw.Balloon(self.chosen_own_link, relmouse="both")
            self.create_type_loop_closing_interior_hints(hint_choose_own_link)

    def create_type_loop_closing_interior_hints(self, current_page):
        current_page.bind(self.btns_get_data.button(0), textwrap.fill("Retrieve the selected tuple of C-alpha atoms "
                                                                      "from the (bio)polymer loaded in PyMOL.",
                                                                      self.hint_width))
        current_page.bind(self.btns_get_data.button(1), textwrap.fill("Retrieve the selected sequence of C-alpha atoms "
                                                                      "from the (bio)polymer loaded in PyMOL.",
                                                                      self.hint_width))
        current_page.bind(self.btns_bottom.button(0), "Clear the above data.")

        if hasattr(self, "btn_whole_chains") and self.btn_whole_chains.winfo_exists():
            current_page.bind(self.btn_whole_chains, textwrap.fill("Calculate the linking between whole chains. "
                                                                   "Alternativelly, one can determine the linking "
                                                                   "between parts of chain specifying the fragment "
                                                                   "endpoints in table below.", self.hint_width))

        for i in self.btns_view:
            current_page.bind(i, "View the loop closed by the the bridge chosen.")

        for i in self.loops_list:
            current_page.bind(i[0], textwrap.fill("Specify chain and the residues delimiting fragment of interest. "
                                                  "Alternativelly tick the option ''Use whole chains'' to investigate "
                                                  "linking between whole chains.", self.hint_width))

    def modify_whole_chains(self):
        if not self.whole_chains_enabled.get():
            for idx, pos in enumerate(self.loops_list):
                pos[1].configure(entry_state='normal')
                pos[2].configure(entry_state='normal')
                self.btns_view[idx].configure(state='normal')
        else:
            for idx, pos in enumerate(self.loops_list):
                pos[1].configure(entry_state='disabled')
                pos[2].configure(entry_state='disabled')
                pos[1].setentry("")
                pos[2].setentry("")
                self.btns_view[idx].configure(state='disabled')

    def clear_loops_list(self):
        for i in self.loops_list:
            i[0].component("entryfield").setentry("")
            i[1].setentry("")
            i[2].setentry("")

        for i in self.distance_error.keys():
            self.distance_error[i].grid_forget()

    def pymol_display_chain(self, list_chains):
        chains = []

        if list_chains:
            for chain in list_chains:
                if ".xyz" in chain:  # ions
                    file = chain.split(".pdb_")[0]
                    chain = chain.split(".pdb_")[1].split("_")[0]
                    chains.append(file + " " + chain)
                else:
                    chain = chain.split(" ")
                    if len(chain) == 3:
                        chains.append(chain[0][:-1] + " " + chain[0][-1])
                    else:
                        chains.append(chain[0] + " " + chain[1])

            cmd.hide(representation="everything", selection="all")
            for pyobj in cmd.get_names(type="selections"):
                if pyobj not in chains:
                    cmd.delete(name=pyobj)

            for chain in chains:
                link_file, chain = chain.split(" ")

                cmd.select(name=link_file + "_" + chain, selection=link_file + " and chain " + chain)
                if self._is_file_xyz[link_file]:
                    cmd.show(representation="lines", selection=link_file + "_" + chain)
                    if self.previously_displayed_bond:
                        if abs(int(self.previously_displayed_bond[1]) - int(self.previously_displayed_bond[2])) >= 2:
                            cmd.unbond(atom1="chain " + self.previously_displayed_bond[0] + " and residue " +
                                             self.previously_displayed_bond[1] + " and name ca",
                                       atom2="chain " + self.previously_displayed_bond[0] + " and residue " +
                                             self.previously_displayed_bond[2] + " and name ca")
                else:
                    cmd.show(representation="cartoon", selection=link_file + "_" + chain)
                cmd.cartoon(type="tube", selection=link_file + "_" + chain)
                util.cbc()
                cmd.orient(selection=link_file + "_*")
                cmd.center(selection=link_file + "_*")
            cmd.deselect()

    def view_possible_loop(self, loop, idx):
        if any([len(_.getvalue()) == 0 for _ in loop]):
            self.raise_popup_menu('There are missing data in the fields.')

        link_file, chain = loop[0].getvalue()[0].split(" ")
        loop_beg = loop[1].getvalue()
        loop_end = loop[2].getvalue()
        atom = "ca"
        if link_file + "_" + chain not in cmd.get_names():
            self.pymol_display_chain([link_file + " " + chain])

        atoms = {'atoms': []}
        cmd.iterate_state(state=0, selection=link_file + " and chain " + chain, expression="atoms.append(resi)",
                          space=atoms)
        atoms = atoms['atoms']
        marginal_end = min(self._chains_marginal_atoms[link_file + chain][1], atoms[-1])
        if (int(loop_beg) < int(atoms[0])) or (int(loop_end) > int(marginal_end)):
            self.raise_popup_menu('There are no such indices in the structure. Please choose other (correct) ones.')

        for i in cmd.get_names(type="all"):
            if str(i).startswith("TMP") or str(i).startswith("SEQ") or str(i).startswith("TRIANG") or \
                    str(i).startswith("BR*") or str(i).startswith("SMOOTH_CHAIN_") or str(i).startswith("SEQ_SM") or \
                    str(i).startswith("SEQ_SM_") or str(i).startswith(link_file + "_*"):
                cmd.hide(representation="everything", selection=i)
                cmd.delete(name=i)

        if self._is_file_xyz[link_file] or not self._contains_bridge[link_file]:
            if self.previously_displayed_bond:
                if abs(int(self.previously_displayed_bond[1]) - int(self.previously_displayed_bond[2])) >= 2:
                    cmd.unbond(atom1="chain " + self.previously_displayed_bond[0] + " and residue " +
                                     self.previously_displayed_bond[1] + " and name " + atom,
                               atom2="chain " + self.previously_displayed_bond[0] + " and residue " +
                                     self.previously_displayed_bond[2] + " and name " + atom)
            self.previously_displayed_bond = [chain, loop_beg, loop_end]
            cmd.show(representation="lines", selection=link_file)
            cmd.set(name="line_width", value="4")

        cmd.set(name="seq_view", value=1)
        cmd.set(name="mouse_selection_mode", value=6)
        cmd.hide(representation="sticks", selection="all")
        cmd.show(representation="cartoon", selection=link_file + "_*")
        cmd.cartoon(type="tube", selection=link_file + "_*")
        util.cbc()
        cmd.select("TMP_BR_" + loop_beg + "_" + loop_end,
                   selection="(chain " + chain + " and residue " + loop_beg + " and name " + atom + ")+" +
                             "(chain " + chain + " and residue " + loop_end + " and name " + atom + ")")
        cmd.select("TMP_SEQ",
                   selection="chain " + chain + " and residue " + loop_beg + "-" + loop_end)
        cmd.bond(atom1="chain " + chain + " and residue " + loop_beg + " and name " + atom,
                 atom2="chain " + chain + " and residue " + loop_end + " and name " + atom)
        cmd.color(color="gray", selection="TMP_SEQ")
        cmd.show(representation='sphere', selection="TMP_BR_" + loop_beg + "_" + loop_end)
        cmd.show(representation="sticks", selection="TMP_BR_" + loop_beg + "_" + loop_end)
        cmd.center(selection="TMP_SEQ*")
        cmd.deselect()

        self.check_distance_between_atoms(loop_beg, loop_end, chain, idx)

    def check_distance_between_atoms(self, beg, end, chain, pos):
        atom = "ca"

        selection1 = "(chain " + chain + " and residue " + beg + " and name " + atom + ")"
        selection2 = "(chain " + chain + " and residue " + end + " and name " + atom + ")"

        cmd.distance(name="DIST_" + beg + "_" + end, selection1=selection1, selection2=selection2)

        calc_distance = round(cmd.get_distance(atom1=selection1, atom2=selection2), 1)
        print "  ## DISTANCE between atoms " + beg + " and " + end + " is " + str(calc_distance)

        if float(calc_distance) > 10:
            if pos not in self.distance_error.keys():
                self.distance_error[pos] = (tk.Label(self.chosen_own_link, text='too big!', foreground="red"))
                self.distance_error[pos].grid(column=5, row=pos + 3)
                self.list_hint_distance.append(Pmw.Balloon(self.chosen_own_link, relmouse="both"))
                self.list_hint_distance[-1].bind(self.distance_error[pos],
                                                 textwrap.fill("A distance between residues is too big to form a "
                                                               "cysteine bridge. Another pair of residues forming "
                                                               "a bridge should be chosen or an option 'Ignore "
                                                               "bad length of bridge or Ca-Ca bonds' should be "
                                                               "ticked.", self.hint_width))
        else:
            if pos in self.distance_error.keys():
                self.distance_error[pos].grid_forget()
                del self.distance_error[pos]

    def get_loop_from_protein_structure(self):
        if "sele" not in cmd.get_names(type="selections"):
            self.raise_popup_menu('No atoms have been chosen.')

        chain = cmd.get_chains(selection="sele")[0]
        link_file = cmd.get_names(selection="sele")[0]

        atoms = {'atoms': []}
        cmd.iterate_state(state=0, selection="sele and chain " + chain + " and name ca",
                          expression='atoms.append(resi)', space=atoms)
        atoms = {'atoms': list(set(atoms['atoms']))}
        atoms['atoms'].sort(key=int)

        if len(atoms['atoms']) > 2:
            self.raise_popup_menu('Too many selection (' + str(len(atoms['atoms'])) + ').')
        elif len(atoms['atoms']) < 2:
            self.raise_popup_menu('Too little selections or there are selected atoms from wrong chain.')

        for idx, loop in enumerate(self.loops_list):
            if len(loop[0].get()) >= 3 and loop[0].get().split(" ")[1] == chain and (not loop[1].get() and
                                                                                         not loop[2].get()):
                self.btn_whole_chains.deselect()
                self.modify_whole_chains()
                loop[1].setentry(atoms['atoms'][0])
                loop[2].setentry(atoms['atoms'][1])
                break
            if not loop[0].get() and not loop[1].get() and not loop[2].get():
                self.btn_whole_chains.deselect()
                self.modify_whole_chains()
                loop[0].selectitem(link_file + " " + chain)
                loop[1].setentry(atoms['atoms'][0])
                loop[2].setentry(atoms['atoms'][1])
                break
        else:
            self.raise_popup_menu('All the fields are filled. There is no space to enter another loop.')

        if len(cmd.get_names(type="selections")) > 0:
            cmd.delete(name="sele")

    def get_loop_from_protein_sequence(self):
        if "sele" not in cmd.get_names(type="selections"):
            self.raise_popup_menu('No atoms have been chosen.')

        chain = cmd.get_chains(selection="sele")[0]
        link_file = cmd.get_names(selection="sele")[0]

        atoms = {'atoms': []}
        cmd.iterate_state(state=0, selection="sele and chain " + chain + " and name ca",
                          expression='atoms.append(resi)', space=atoms)
        atoms = {'atoms': list(set(atoms['atoms']))}
        atoms['atoms'].sort(key=int)

        if len(atoms['atoms']) < 2:
            self.raise_popup_menu('Too little selections or there are selected atoms from wrong chain.')

        for idx, loop in enumerate(self.loops_list):
            if len(loop[0].get()) >= 3 and loop[0].get().split(" ")[1] == chain and (not loop[1].get() and
                                                                                         not loop[2].get()):
                self.btn_whole_chains.deselect()
                self.modify_whole_chains()
                loop[1].setentry(atoms['atoms'][0])
                loop[2].setentry(atoms['atoms'][1])
                break
            if not loop[0].get() and not loop[1].get() and not loop[2].get():
                self.btn_whole_chains.deselect()
                self.modify_whole_chains()
                loop[0].selectitem(link_file + " " + chain)
                loop[1].setentry(atoms['atoms'][0])
                loop[2].setentry(atoms['atoms'][1])
                break
        else:
            self.raise_popup_menu('All the fields are filled. There is no space to enter another loop.')

        if len(cmd.get_names(type="selections")) > 0:
            cmd.delete(name="sele")

    def choose_ions_closing(self):
        if any(self._is_file_xyz[elem] for elem in self._is_file_xyz):
            self.choose_automatic_closing()
            self.raise_popup_menu("XYZ files contain no information about chemical bridges, therefore "
                                  "calculation of extended-type loops is not available.")

        if hasattr(self, "chosen_own_link"):
            self.chosen_own_link.grid_forget()
            self.type_closing_choice[0] = False

        if not self.type_closing_choice[1]:
            self.type_closing_choice[1] = True
            self.is_ions = True
            self.chosen_own_ions_link = Pmw.ScrolledFrame(self.group_method.interior(), hull_width=500, hull_height=400,
                                                          usehullsize=1)
            self.chosen_own_ions_link.grid(sticky='swen', column=0, row=2, padx=5, pady=5)

            extended_label = tk.Label(self.chosen_own_ions_link.interior(),
                                      text=textwrap.fill('Pick ion loops from which you want to create a link. Up to '
                                                         'four loops are available.', 80))
            extended_label.grid(column=0, row=0, columnspan=5, pady=2, padx=5)

            if not hasattr(self, "ions_connections"):
                self.retrieve_ions_data()
            if len(self.ions_connections) >= 2:
                self.create_ions_links_buttons()
            else:
                loops_list_label = tk.Label(self.chosen_own_ions_link.interior(),
                                            text=textwrap.fill('The number of links found in given structure is not '
                                                               'enough to form any link.', 70))
                loops_list_label.grid(column=0, row=2, columnspan=5, pady=2, padx=5)
                self.chosen_own_ions_link.configure(hull_height=500)

    def retrieve_ions_data(self):
        self.ions_connections = {}
        args = []

        for filename in self._full_path_to_files:
            file = filename.split(os.sep)[-1][:-4]
            filtered_ions = filter(len, os.popen(self.prog_converter + filename + " -e").read().splitlines())
            re_atoms = re.compile("\w*-\d_\w\d*")

            if filtered_ions and ("ERROR!!!" in filtered_ions[0] or "WARNING!!!" in filtered_ions[0]):
                self.raise_popup_menu(filtered_ions[0])

            elif filtered_ions:
                for line in filtered_ions:
                    file_number = line.split("\t")[0]
                    filename = line.split("\t")[1].split(os.sep)[-1]
                    loop = " ".join(line.split("\t")[2:])
                    atoms = re_atoms.findall(loop)
                    cysteine = atoms[:2]
                    first_atom = int(loop.split(" ")[0].split("_")[1][1:])
                    last_atom = int(loop.split(" ")[-3].split("_")[1][1:])
                    relation = loop.split(" ")[-2]
                    bridge = file + cysteine[0].split("_")[1][0] + " "
                    if relation == "<->":
                        bridge += str(first_atom) + " " + str(last_atom)
                    elif last_atom > first_atom:
                        bridge += str(first_atom) + " " + str(first_atom + 1)
                    else:
                        bridge += str(first_atom) + " " + str(first_atom - 1)
                    args.append(file_number)
                    self.ions_connections[int(file_number)] = [filename, bridge, loop]

        return args

    def create_ions_links_buttons(self):
        self.ion_lines = []
        self.is_ion_link_clicked = []

        for idx, link in enumerate(sorted(self.ions_connections.keys())):
            atoms = self.ions_connections[link][2]
            content = tk.Text(self.chosen_own_ions_link.interior(), bg="white", padx=0, pady=0, wrap="word",
                              bd=0, highlightthickness=0, width=72, height=textwrap.fill(atoms, 70).count("\n") + 2)
            self.color_ion_loop_text(content, atoms, 0)
            content.grid(row=idx + 1, column=1)

            self.is_ion_link_clicked.append(tk.IntVar())
            tick = tk.Checkbutton(self.chosen_own_ions_link.interior(), width=2, text=str(idx) + " ", variable=self.is_ion_link_clicked[-1],
                                  command=lambda x=idx: self.select_ion_loop(x))
            tick.grid(row=idx + 1, column=0)
            self.ion_lines.append([tick, content])

    def color_ion_loop_text(self, tk_txt_elem, txt, height):
        tk_txt_elem.configure(state="normal")
        amino_acids = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET',
                       'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL', 'BTC', 'FCY', 'GGL']

        col = 0.0
        beg = 0
        elem_len = 0
        for line in txt.split(" "):
            col += 1.0
            tk_txt_elem.insert(col, line + " ")
            beg += elem_len
            elem_len = len(line) + 1
            if "<->" not in line and "..." not in line:
                if "CYS" in line:
                    color = "orange"
                elif line.split("-")[0] in amino_acids:
                    color = "blue"
                else:
                    color = "black"
                tk_beg = "1." + str(beg)
                tk_end = "1." + str(beg + elem_len)
                tk_txt_elem.tag_add(color, tk_beg, tk_end)
                tk_txt_elem.tag_config(color, foreground=color)
        tk_txt_elem.tag_add("center", 1.0, "end")
        tk_txt_elem.tag_configure("center", justify='center')
        tk_txt_elem.tag_add("spacing1", 1.0, "end")
        tk_txt_elem.tag_config("spacing1", spacing1=5 + height)
        tk_txt_elem.configure(state="disabled")

    def select_ion_loop(self, idx):
        if self.is_ion_link_clicked[idx].get():
            if [i.get() for i in self.is_ion_link_clicked].count(1) > 4:
                self.ion_lines[idx][0].deselect()
                self.raise_popup_menu('Only four components can be chosen. Unclick one to select next one.')

            self.ion_lines[idx][1].configure(bg="lightblue")
            self.prepare_pymol_ion_structures(idx)
        else:
            self.ion_lines[idx][1].configure(bg="white")

    def prepare_pymol_ion_structures(self, idx):
        cmd.hide(representation="everything", selection="all")
        for pyobj in cmd.get_names(type="selections"):
            cmd.delete(name=pyobj)
        util.cbc()
        cmd.hide(representation="spheres", selection="all")
        self.pymol_display_ion_structures(self.ions_connections[idx][2].split(" "))

    def pymol_display_ion_structures(self, atoms, color_ion_selection=True):
        amino_acids = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET',
                       'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL', 'BTC', 'FCY', 'GGL']

        for at1, rel, at2 in zip(atoms[::2], atoms[1::2], atoms[2::2]):
            at1 = at1.split("_")
            at2 = at2.split("_")

            if int(at1[1][1:]) > int(at2[1][1:]): # if the first residue is bigger than second, we have to swap them
                tmp = at1
                at1 = at2
                at2 = tmp

            if at1[0].split("-")[0] in amino_acids:
                first_atom = "chain " + at1[1][0] + " and residue " + at1[1][1:] + " and resn " + \
                             at1[0].split("-")[0] + " and name ca"
                first_color = "orange"
            else:
                first_atom = "chain " + at1[1][0] + " and residue " + at1[1][1:] + " and resn " + \
                             at1[0].split("-")[0]
                first_color = "blue"
            if at2[0].split("-")[0] in amino_acids:
                second_atom = "chain " + at2[1][0] + " and residue " + at2[1][1:] + " and resn " + \
                              at2[0].split("-")[0] + " and name ca"
                second_color = "orange"
            else:
                second_atom = "chain " + at2[1][0] + " and residue " + at2[1][1:] + " and resn " + \
                              at2[0].split("-")[0]
                second_color = "blue"
            marginal_atoms = "(" + first_atom + ")+(" + second_atom + ")"

            if rel == "...":
                cmd.show(representation="cartoon", selection="chain " + at1[1][0])
                cmd.show(representation='sphere', selection=marginal_atoms)
                if color_ion_selection:
                    cmd.color(color="radium", selection="chain " + at1[1][0] + " and residue " + at1[1][1:] +
                                                        "-" + at2[1][1:] + " and name ca")
                cmd.show(representation='sticks', selection=marginal_atoms)
            else:
                cmd.show(representation='sphere', selection=marginal_atoms)
                cmd.show(representation="sticks", selection=marginal_atoms)
                cmd.bond(atom1=first_atom, atom2=second_atom)
                cmd.set(name="stick_color", value="blue", selection="all")

                cmd.set(name="sphere_color", value=first_color, selection=first_atom)
                cmd.set(name="sphere_color", value=second_color, selection=second_atom)
        cmd.center(selection=cmd.get_names()[0])
        cmd.deselect()

    def create_smooth_interior(self):
        self.group_smooth = Pmw.Group(self.group_advanced.interior(), tag_text='Smoothing')
        self.group_smooth.grid(sticky='eswn', column=1, row=0, padx=5, pady=5)

        self.smooth_level = Pmw.EntryField(self.group_smooth.interior(), labelpos='w', label_text='Level of smoothness',
                                           validate={'validator': 'integer', 'min': 1, "max": 101},
                                           value=2, entry_width=6)
        self.smooth_level.grid(sticky='w', column=0, row=0, padx=10, pady=2)

    def create_structure_correctness_interior(self):
        self.group_structure_correctness = Pmw.Group(self.group_advanced.interior(), tag_text='Correct mode')
        self.group_structure_correctness.grid(sticky='eswn', column=1, row=1, padx=5, pady=5)

        self.structure_correctness = tk.Checkbutton(self.group_structure_correctness.interior(),
                                                    text='Check structure correctness', variable=self.is_correct_data)
        self.structure_correctness.grid(sticky='w', column=0, row=1, padx=2, pady=2)

    def create_gln_interior(self):
        self.group_gln = Pmw.Group(self.group_advanced.interior(), tag_text='Gaussian Linking Number')
        self.group_gln.grid(sticky='eswn', column=1, row=3, padx=5, pady=5)

        gln_checkbutton = tk.Checkbutton(self.group_gln.interior(), text='Calculate GLN', variable=self.is_gln_selected)
        gln_checkbutton.grid(sticky='w', column=0, row=0, padx=2, pady=2)

    def create_advanced_frame_hints(self, current_page):
        current_page.bind(self.group_smooth, textwrap.fill("Adjust the level of smoothness, by choosing the value of "
                                                           "this parameter in the range 2-100. The higher value "
                                                           "results in more refined triangulation of the minimal "
                                                           "surface. Above a certain value (characteristic for a "
                                                           "given configuration) there is no visible change in the "
                                                           "triangulation.", self.hint_width))
        current_page.bind(self.group_structure_correctness, textwrap.fill("Tick this option to perform basic "
                                                                          "structure validation. Structures with "
                                                                          "unnatural distance between C-alpha atoms "
                                                                          "will be rejected.", self.hint_width))
        current_page.bind(self.group_gln, textwrap.fill("A numerical invariant that describes the linking of two "
                                                        "closed curves in the three-dimensional space. "
                                                        "Intuitively, the linking number represents the number of "
                                                        "times that each curve winds around the other one.",
                                                        self.hint_width))
        if hasattr(self, "group_number_of_projections"):
            current_page.bind(self.group_number_of_projections, textwrap.fill("The number of different chain closures "
                                                                              "to calculate. Higher number enhances the"
                                                                              " precision of the link probability "
                                                                              "calculation.", self.hint_width))

    def create_probabilistic_page(self):
        if hasattr(self, "disable_pages") and self.disable_pages:
            self.notebook.selectpage("Macromolecular links")
            self.raise_popup_menu("No probabilistic approach for macromolecular links.")

        self.hide_interface_components()

        self.page_probabilistic = tk.Label(self.notebook.page(1))
        self.page_probabilistic.grid(sticky='eswn', column=0, row=0)

        self.group_method = Pmw.Group(self.page_probabilistic, tag_text='Link identification method')
        self.group_method.grid(sticky='eswn', column=0, row=0, padx=5, pady=5)

        self.identification_method = Pmw.ComboBox(self.group_method.interior(), entry_width=30,
                                                  scrolledlist_items=('automatic detection',
                                                                      'choose bridging atoms'),
                                                  entryfield_entry_state="readonly",
                                                  selectioncommand=lambda x: self.set_link_closing())
        self.identification_method.selectitem(0)
        self.identification_method.grid(sticky='eswn', column=0, row=1, padx=5, pady=5)

        self.group_advanced = Pmw.Group(self.page_probabilistic, tag_text='Advanced options')
        self.group_advanced.grid(sticky='eswn', column=1, row=0, padx=5, pady=5, rowspan=3)

        self.create_smooth_interior()
        self.create_structure_correctness_interior()
        self.create_gln_interior()
        self.create_number_of_projections_interior()

        hint_advanced = Pmw.Balloon(self.page_probabilistic, relmouse="both")
        self.create_advanced_frame_hints(hint_advanced)

        self.notebook.setnaturalsize()

    def create_number_of_projections_interior(self):
        self.group_number_of_projections = Pmw.Group(self.group_advanced.interior(), tag_text='Number of tries')
        self.group_number_of_projections.grid(sticky='eswn', column=1, row=4, padx=5, pady=5)

        self.number_of_projections = Pmw.EntryField(self.group_number_of_projections.interior(), labelpos='w',
                                                    label_text='Number of tries', validate={'validator': 'integer',
                                                                                            'min': 1, "max": 1001},
                                                    value=100, entry_width=5)
        self.number_of_projections.grid(sticky='w', column=1, row=1, padx=10, pady=2)

    def create_macromolecular_page(self):
        self.hide_interface_components()

        self.page_macromolecular = tk.Label(self.notebook.page(2))
        self.page_macromolecular.grid(sticky='eswn', column=0, row=0)

        self.group_method = Pmw.Group(self.page_macromolecular, tag_text='Link identification method')
        self.group_method.grid(sticky='eswn', column=0, row=0, padx=5, pady=5)

        self.identification_method = Pmw.ComboBox(self.group_method.interior(), entry_width=30,
                                                  scrolledlist_items=('automatically join chains',
                                                                      'choose joining residues'),
                                                  entryfield_entry_state="readonly",
                                                  selectioncommand=lambda x: self.set_macrolink_closing())
        self.identification_method.grid(sticky='eswn', column=0, row=1, padx=5, pady=5)
        self.identification_method.selectitem(0)
        self.set_macrolink_closing()

        self.group_advanced = Pmw.Group(self.page_macromolecular, tag_text='Advanced options')
        self.group_advanced.grid(sticky='eswn', column=1, row=0, padx=5, pady=5, rowspan=3)

        self.create_smooth_interior()
        self.create_structure_correctness_interior()
        self.create_gln_interior()

        hint_advanced = Pmw.Balloon(self.page_macromolecular, relmouse="both")
        self.create_advanced_frame_hints(hint_advanced)

        self.notebook.setnaturalsize()

    def set_macrolink_closing(self):
        if str(self.identification_method.get()) == 'automatically join chains':
            self.choose_automatically_join_chains()
        else:
            self.choose_joining_residues()
        self.notebook.setnaturalsize()

    def choose_automatically_join_chains(self):
        if hasattr(self, "chosen_own_link"):
            self.chosen_own_link.grid_forget()
            self.type_closing_choice[1] = False

        if not self.type_closing_choice[0]:
            self.type_closing_choice[0] = True
            cmd.set(name="mouse_selection_mode", value=2)

            self.chosen_automatic_macrolink = tk.Frame(self.group_method.interior())
            self.chosen_automatic_macrolink.grid(sticky='swen', column=0, row=2, padx=5, pady=5)

            macrolink_label = tk.Label(self.chosen_automatic_macrolink,
                                       text=textwrap.fill("Specify chains forming each component. Each part must "
                                                          "consist of model's number and chain binded with '_' sign "
                                                          "(e.g. 1_A, 2_B, 3_C, 4_D). Remember about spacing.", 50))
            macrolink_label.grid(column=0, row=0, columnspan=4, pady=2, padx=5)

            self.macrolink_list = []
            component_labels = {0: "First component     ", 1: "Second component", 2: "Third component    ",
                                3: "Fourth component  "}
            for i in range(4):
                self.macrolink_list.append(Pmw.EntryField(self.chosen_automatic_macrolink, labelpos='w',
                                                          label_text=component_labels[i], entry_width=25,
                                                          validate=None))
                self.macrolink_list[-1].grid(column=0, row=i + 3, pady=10)

            self.btns_bottom = Pmw.ButtonBox(self.chosen_automatic_macrolink, orient="horizontal")
            self.btns_bottom.add("Clear", width=self.bottom_btns_width, command=lambda: self.clear_macrolinks_list())
            self.btns_bottom.add("Import", width=self.bottom_btns_width, command=lambda: self.import_macrocomponents())
            self.btns_bottom.grid(sticky='swn', column=0, row=len(self.macrolink_list) + 3, columnspan=2, padx=5, pady=2)

            self.create_joining_method_interior()

            hint_choose_automatic_macrolinks = Pmw.Balloon(self.chosen_automatic_macrolink, relmouse="both")
            self.create_macrolink_automatic_hints(hint_choose_automatic_macrolinks)
        self.notebook.setnaturalsize()

    def clear_macrolinks_list(self):
        for i in self.macrolink_list:
            i.setentry("")

    def import_macrocomponents(self,):
        self.macrocomponent_file = tkFileDialog.askopenfilenames(initialdir=os.getcwd(), filetypes=[("TXT", "*.txt")],
                                                                 title="Choose file with macrocomponents")

        if not self.macrocomponent_file or None == self.macrocomponent_file:
            return
        elif len(self.macrocomponent_file) > 1:
            print "  ### Too many files were chosen."
            return
        else:
            try:
                content = open(self.macrocomponent_file[0], "r")
                for i, line in enumerate(content):
                    line = line.rstrip()
                    if line != "\n":
                        self.macrolink_list[i].setentry(line)
            except IndexError:
                self.raise_popup_menu("An error occurred during import of file. Please make sure the lines contain "
                                      "proper data.")
        self.notebook.setnaturalsize()

    def create_joining_method_interior(self):
        self.group_joining_method = Pmw.Group(self.chosen_automatic_macrolink, tag_text='Method to join chains')
        self.group_joining_method.grid(sticky='eswn', column=0, row=len(self.macrolink_list) + 4, columnspan=2, padx=5,
                                       pady=5)

        self.method_to_join_chains = Pmw.RadioSelect(self.group_joining_method.interior(), buttontype='radiobutton',
                                                     orient='vertical')
        self.method_to_join_chains.grid(sticky='w', column=0, row=0, padx=2, pady=2)
        for idx, label in enumerate(["Join nearest residues", "Join termini"]):
            self.method_to_join_chains.add(label)
        self.method_to_join_chains.invoke("Join nearest residues")

    def create_macrolink_automatic_hints(self, current_page):
        current_page.bind(self.btns_bottom.button(0), "Clear the above data.")
        current_page.bind(self.btns_bottom.button(1), "Import residues creating component from file.txt.")
        current_page.bind(self.method_to_join_chains, "Select method to join chains either by nearest residue or termini.")

    def choose_joining_residues(self):
        if hasattr(self, "chosen_automatic_macrolink"):
            self.chosen_automatic_macrolink.grid_forget()
            self.type_closing_choice[0] = False

        if not self.type_closing_choice[1]:
            self.macrolink_list = []
            self.component_groups = []
            self.btns_bottom = []
            self.type_closing_choice[1] = True
            cmd.set(name="mouse_selection_mode", value=6)

            self.chosen_own_link = tk.Frame(self.group_method.interior())
            self.chosen_own_link.grid(sticky='swen', column=0, row=2, padx=5, pady=5)

            loops_list_label = tk.Label(self.chosen_own_link, text=textwrap.fill('Specify the bridges between chains '
                                                                                 'creating a macrocomponent. Each atom'
                                                                                 ' forming a bridge must contain '
                                                                                 'information about the model, the '
                                                                                 'chain, the index of the residue and '
                                                                                 'be separated with a "_" character '
                                                                                 '(e.g. 1_A_56).', 70))
            loops_list_label.grid(column=0, row=0, columnspan=5, pady=2, padx=5)

            for idx, comp in enumerate(["First component", "Second component", "Third component", "Fourth component"]):
                self.component_groups.append(Pmw.Group(self.chosen_own_link, tag_text=comp))
                if idx > 1:
                    self.component_groups[-1].grid(sticky="swen", column=idx-2, row=2, padx=2, pady=2)
                else:
                    self.component_groups[-1].grid(sticky="swen", column=idx, row=1, padx=2, pady=2)
                self.create_own_component_interior(self.component_groups[-1].interior(), idx)

            self.create_joining_residues_hints()
            self.notebook.setnaturalsize()

    def create_own_component_interior(self, cmpt_interior, idx):
        comp_list = []

        for j in range(5):
            comp_list.append((
                (Pmw.EntryField(cmpt_interior, entry_width=7, validate=None)),
                (Pmw.EntryField(cmpt_interior, entry_width=7, validate=None))))
            comp_list[-1][0].grid(column=2, row=j + 1, padx=5, pady=5)
            comp_list[-1][1].grid(column=3, row=j + 1, padx=5, pady=5)

        btns = Pmw.ButtonBox(cmpt_interior, orient="vertical")
        btns.add("Extend", width=self.bottom_btns_width, command=lambda x=j: self.extend_macrolink_list(idx))
        btns.add("Import", width=self.bottom_btns_width, command=lambda x=j: self.import_macrolink_list(self.macrolink_list, idx))
        btns.add("Clear", width=self.bottom_btns_width, command=lambda x=j: self.clear_macrolink_list(idx))
        btns.grid(sticky='swn', column=0, row=0, rowspan=5, padx=5, pady=2)
        self.btns_bottom.append(btns)
        self.macrolink_list.append(comp_list)

    def extend_macrolink_list(self, idx):
        macrolist = self.macrolink_list[idx]
        comp_interior = self.component_groups[idx].interior()

        if len(macrolist) == 13:
            self.raise_popup_menu("Only 13 loops can be included in a macrocomponent.")
        if not all(len(elem[0].get()) != 0 and len(elem[1].get()) != 0 for elem in macrolist):
            self.raise_popup_menu("All fields must be filled in order to add a new pair.")

        macrolist.append((Pmw.EntryField(comp_interior, entry_width=7),
                          Pmw.EntryField(comp_interior, entry_width=7)))
        macrolist[-1][0].grid(column=2, row=len(macrolist) + 3, padx=5, pady=5)
        macrolist[-1][1].grid(column=3, row=len(macrolist) + 3, padx=5, pady=5)
        self.notebook.setnaturalsize()

    def import_macrolink_list(self, macrolist, idx):
        self.macrocomponent_file = tkFileDialog.askopenfilenames(initialdir=os.getcwd(), filetypes=[("TXT", "*.txt")],
                                                                 title="Choose a file with macrocomponent")

        if not self.macrocomponent_file or self.macrocomponent_file == None:
            return
        elif len(self.macrocomponent_file) > 1:
            print "  ### Too many files were chosen."
            return
        else:
            content = open(self.macrocomponent_file[0], "r")
            for i, line in enumerate(content):
                line = line.rstrip().split(" ")
                if len(line) != 1:
                    if i >= len(self.macrolink_list[idx]):
                        macrolist[idx].append((Pmw.EntryField(self.component_groups[idx].interior(), entry_width=7),
                                               Pmw.EntryField(self.component_groups[idx].interior(), entry_width=7)))
                        macrolist[idx][-1][0].grid(column=2, row=len(macrolist[idx]) + 3, padx=5, pady=5)
                        macrolist[idx][-1][1].grid(column=3, row=len(macrolist[idx]) + 3, padx=5, pady=5)
                    macrolist[idx][i][0].setentry(line[0])
                    macrolist[idx][i][1].setentry(line[1])
        self.notebook.setnaturalsize()

    def clear_macrolink_list(self, idx):
        for i in self.macrolink_list[idx][:5]:
            i[0].setentry("")
            i[1].setentry("")
        for i, elem in enumerate(self.macrolink_list[idx][5:]):
            elem[0].grid_forget()
            elem[1].grid_forget()
        self.macrolink_list[idx] = self.macrolink_list[idx][:5]
        self.notebook.setnaturalsize()

    def create_joining_residues_hints(self):
        hint_joining_residues = Pmw.Balloon(self.group_joining_method.interior(), relmouse="both")

        for elem in self.btns_bottom:
            hint_joining_residues.bind(elem.button(0), textwrap.fill("Extend the current list.", self.hint_width))
            hint_joining_residues.bind(elem.button(1), textwrap.fill("Import residues creating this component from "
                                                                     "file.txt."
                                                                     , self.hint_width))
            hint_joining_residues.bind(elem.button(2), textwrap.fill("Clear the list.", self.hint_width))

    def _invoke_plugin_action(self, clicked_btn):
        if clicked_btn == "Proceed":
            print "  PyLink is running..."
            self._invoke_program()
        else:
            self.dialog.withdraw()
            if hasattr(self, "error_pop_menu") and self.error_popup.winfo_exists():
                self.error_popup.withdraw()
            if hasattr(self, "win_link_info") and self.win_link_info.winfo_exists():
                self.win_link_info.withdraw()
            if hasattr(self, "artifact_found") and self.artifact_found.winfo_exists():
                self.artifact_found.withdraw()
            if hasattr(self, "ions_connections") and self.ions_connections and hasattr(self, "link_directory"):
                for filename in self._filenames:
                    self.separate_files_to_directory(self._full_path_to_dir, self.link_directory,
                                                     filename + "_\w*_\d*.(xyz|pdb)")
            print "  PyLink has been shut down..."

    ####################################################################################################################
    #                                               EXECUTE PROGRAM 1/2
    ####################################################################################################################

    def _invoke_program(self):
        self.displayed_link = None
        self.output_data = {}

        self.convert_to_5columns_format()
        if self.is_correct_data.get() and not all(len(elem) == 0 for elem in self.warning_gaps):
            self.raise_popup_menu("Given structure has gaps or there are unnatural distances between atoms. "
                                  "Calculations aborted.")
        if self.identification_method.get() == "choose bridging atoms (extended on ions)" and \
                len([idx for idx, elem in enumerate(self.is_ion_link_clicked) if elem.get()]) < 2:
            self.raise_popup_menu("At least two loops are required to form a link.")

        self.purge_pymol_structures()
        self.user_data = self.generate_invoking_commands()
        self.create_main_link_directory()
        self.find_links()
        self.get_links_intersections()
        self.separate_link_files_to_hashdirectories()
        self.purge_output_data()
        self.check_smoothing_avaliability()
        self.purge_main_directory()
        self.calculate_smoothed_structures()

        for filename in self._filenames:
            self.separate_files_to_directory(self._full_path_to_dir, self.link_directory,
                                             filename + "(_[\d|\w].(xyz|pdb))|(macro_(\d\w(_)*)*)")

        if hasattr(self, "win_link_info") and self.win_link_info.winfo_exists():
            self.win_link_info.destroy()
        self.create_link_window()
        if hasattr(self, "error_pop_menu") and self.error_popup.winfo_exists():
            self.error_popup.withdraw()

    def convert_to_5columns_format(self):
        self.pdb_bridges = []
        self.warning_gaps = []

        for filename in self._full_path_to_files:
            file = filename.split(os.sep)[-1][:-4]
            filtered_bridges = filter(len, os.popen(self.prog_converter + filename + " -f").read().splitlines())
            if not self._is_file_xyz[file]:
                self.warning_gaps.append([(i.split(" ")[3], i.split(" ")[-3], i.split(" ")[-1]) for i in
                                          filtered_bridges if "WARNING" in i])
                self.pdb_bridges.append([num for num in filtered_bridges if "WARNING" not in num])
            else:
                self.warning_gaps.append([])
                self.pdb_bridges.append([])

    def purge_pymol_structures(self):
        for pyobj in cmd.get_names(type="all"):
            if pyobj.startswith("SURF_") or pyobj.startswith("PIERC_") or pyobj.startswith("BR_") \
                    or pyobj.startswith("SM") or pyobj.startswith("NEG_") or pyobj.startswith("POS_") or \
                    pyobj.startswith("SEQ") or pyobj.startswith("macro_"):
                cmd.delete(name=pyobj)
        cmd.hide(representation="sticks", selection="all")
        util.cbc()

    ####################################################################################################################
    #                                               1) GET PLUGIN DATA
    ####################################################################################################################

    def generate_invoking_commands(self):
        print "  Collecting user's data..."
        arguments = []

        curr_page = self.notebook.getcurselection()
        if curr_page == "Deterministic":
            arguments += self.get_deterministic_data()
        elif curr_page == "Probabilistic":
            arguments += self.get_probabilistic_data()
        else:
            arguments += self.get_macromolecular_data()
        if len(self.smooth_level.getvalue()) > 0:
            arguments = [i + " -sm" for i in arguments]
            arguments = [i + " " + str(self.smooth_level.getvalue()) for i in arguments]
        if self.is_gln_selected.get():
            arguments = [i + " -glns" for i in arguments]
            arguments = [i + " 1" for i in arguments]
        arguments = [i + " -sf 2" for i in arguments]
        arguments = [i + " -inter 1" for i in arguments]
        if hasattr(self, "ions_connections") and self.ions_connections:
            arguments = [i + " -noCheckIds 1" for i in arguments]

        return arguments

    def get_deterministic_data(self):
        global global_combinations

        if self.identification_method.get() == "automatic detection":
            if any(self._is_file_xyz[elem] for elem in self._is_file_xyz):
                self.raise_popup_menu("No automatic detection of links for .xyz files.")

            chain_tuple = self.get_data_from_pdb_bridges()
            self.generate_chain_permutation(chain_tuple)
            args = self.generate_data("0")
        elif str(self.identification_method.get()) == "choose bridging atoms":
            self.validate_own_loops()
            chain_tuple = self.get_data_from_loops_list()
            self.combinations = [tuple(chain_tuple)]
            global_combinations = self.combinations
            args = self.generate_data("0")
        else:
            chain_tuple = self.get_data_from_ions_list()
            self.combinations = [tuple(chain_tuple)]
            global_combinations = self.combinations
            args = self.generate_ion_data()

        return args

    def get_data_from_pdb_bridges(self):
        chain_tuple = []

        for file_bridges in self.pdb_bridges:
            for elem in file_bridges:
                elem = elem.split(" ")
                file = elem[1].split(os.sep)[-1][:-10]
                chain = elem[1].split(os.sep)[-1][-5]
                bridge = elem[2] + " " + elem[3]
                chain_tuple.append(file + chain + " " + bridge)

        return chain_tuple

    def validate_own_loops(self):
        tmp_minimum = [len(i[0].get()) == 0 and len(i[1].get()) == 0 and len(i[2].get()) == 0 for i in self.loops_list]
        if tmp_minimum.count(False) < 2:
            self.raise_popup_menu('At least two links must be given to form a link.')

        tmp_duplicates = [(i[0].get(), i[1].get(), i[2].get()) for i in self.loops_list]
        if any(tmp_duplicates.count(x) > 1 for x in tmp_duplicates if list(x) != ['', '', '']):
            self.raise_popup_menu('No duplicates allowed. They won\'t create a super-duper link.')

        for idx, line in enumerate(self.loops_list):
            if not self.whole_chains_enabled.get():
                if len(line[0].get()) == 0 and (len(line[1].get()) != 0 or len(line[2].get()) != 0):
                    self.raise_popup_menu(
                        'There are missing data in line ' + str(idx + 1) + ". Please provide a chain.")
                elif len(line[0].get()) != 0 and ((len(line[1].get()) == 0 and len(line[2].get()) != 0) or
                                                      (len(line[1].get()) != 0 and len(line[2].get()) == 0)):
                    self.raise_popup_menu('There are missing data in line ' + str(idx + 1) + ".")

    def generate_chain_permutation(self, chains):
        self.combinations = []
        global global_combinations

        for i in range(2, 5):
            for comb in combinations(chains, i): # where chains is given argument (chain_tuple)
                if len(comb) < 5:
                    self.combinations.append(comb)

        global_combinations = self.combinations

    def generate_data(self, closure_method, num_tries=""):
        args = []
        regex = re.compile("macro_(\d\w(_)*)*")

        for perm in self.combinations:
            files = ""
            for comp in perm:
                dict_key, comp_beg, comp_end = comp.split(" ")
                if regex.match(dict_key):
                    filename = dict_key + ".xyz"
                else:
                    filename = dict_key[:-1] + ".pdb_" + dict_key[-1] + ".xyz"
                files += self._full_path_to_dir + os.sep + filename + " " + comp_beg + " " + comp_end + " "
            n_components = str(len(perm))
            args.append(self.prog_homfly + n_components + " " + closure_method + " " + num_tries + "0 " + files[:-1])

        return args

    def get_data_from_loops_list(self):
        chain_tuple = []

        for idx, elem in enumerate(self.loops_list):
            if len(elem[0].get()) >= 3:
                chain = elem[0].get().split(" ")
                data = chain[0] + chain[1]
                if elem[1].get() == "" and elem[2].get() == "":
                    beg = self._chains_marginal_atoms[data][0]
                    end = self._chains_marginal_atoms[data][1]
                else:
                    beg = elem[1].get()
                    end = elem[2].get()
                chain_tuple.append(data + " " + beg + " " + end)

        return chain_tuple

    def get_data_from_ions_list(self):
        chain_tuple = []

        pos = [idx for idx, elem in enumerate(self.is_ion_link_clicked) if elem.get()]
        for idx in pos:
            bridge = " ".join(self.ions_connections[idx][1].split(" ")[1:])
            chain_tuple.append(self.ions_connections[idx][0] + " " + bridge)

        return chain_tuple

    def generate_ion_data(self):
        args = []

        for perm in self.combinations:
            files = ""
            for num in perm:
                files += self._full_path_to_dir + os.sep + num + " "
            n_components = str(len(perm))
            args.append(self.prog_homfly + n_components + " 0 0 " + files[:-1])

        return args

    def get_probabilistic_data(self):
        if len(self.number_of_projections.get()) == 0:
            self.raise_popup_menu("Wrong number of projections. Please insert appropriate value from range <1;1000>.")

        num_tries = self.number_of_projections.get() + " "
        if self.identification_method.get() == "automatic detection":
            chain_tuple = self.generate_probabilistic_data()
            self.generate_chain_permutation(chain_tuple)
        else:
            self.validate_own_loops()
            chain_tuple = self.get_data_from_loops_list()
            global global_combinations
            self.combinations = [tuple(chain_tuple)]
            global_combinations = self.combinations
        args = self.generate_data("2", num_tries)

        return args

    def generate_probabilistic_data(self):
        chain_tuple = []

        for idx, file in enumerate(self._filenames):
            chains = self._chains[idx]
            for chain in chains:
                filechain = file[:-4] + chain
                chain_marginal_atoms = self._chains_marginal_atoms[filechain]
                chain_tuple.append(filechain + " " + chain_marginal_atoms[0] + " " + chain_marginal_atoms[1])

        return chain_tuple

    def get_macromolecular_data(self):
        global global_combinations

        if self.identification_method.get() == "automatically join chains":
            closure = "-m" if self.method_to_join_chains.getvalue() == "Join nearest residues" else "-n"
            chain_tuple = self.generate_components_files(closure)
            self.combinations = [tuple(chain_tuple)]
            global_combinations = self.combinations
            args = self.generate_data("0")
        else:
            components_macronames = self.generate_bondfile()
            self.convert_bondfile_to_components()
            chain_tuple = self.rename_macrocomponent_files(components_macronames)
            self.separate_files_to_directory(system_working_directory, self._full_path_to_dir, "macro_(\d\w(_)*)*")
            self.combinations = [tuple(chain_tuple)]
            global_combinations = self.combinations
            args = self.generate_data("0")

        return args

    def generate_components_files(self, arg):
        print "  Generating component files..."
        if all(len(ch.get()) == 0 for ch in self.macrolink_list):
            self.raise_popup_menu("There are empty fields in the frame fields. Please fill at least two of them.")
        if not all(len(ch.get()) != 0 for ch in self.macrolink_list[:2]):
            self.raise_popup_menu("At least two components must be given.")

        macrolink = self._filenames[0]
        chain_tuple = []
        reg = re.compile("(\d)*_[\w](,)*")

        for idx, comp in enumerate(self.macrolink_list):
            if comp.get() != "":
                if not reg.match(comp.get()):
                    self.raise_popup_menu("Component " + str(idx + 1) + "contains invalid characters.")

                data = ",".join(sorted(comp.get().replace(",", "").split(" ")))
                command = self.prog_converter + self._full_path_to_dir + os.sep + macrolink + " " + arg + " " + data
                print "Creating a component..."
                res = os.popen(command).read()
                if "ERROR!!!" in res:
                    err_msg = res.split("ERROR!!!")[1]
                    self.raise_popup_menu(err_msg)

                self.separate_files_to_directory(system_working_directory, self._full_path_to_dir, "macro_(\d\w(_)*)*")
                macrofile_name = "macro_" + data.replace("_", "").replace(",", "_") + ".xyz"
                macrofile_atoms = self.get_macrocycle_marginal_atoms(self._full_path_to_dir + os.sep + macrofile_name)
                chain_tuple.append(macrofile_name[:-4] + " " + macrofile_atoms[0] + " " + macrofile_atoms[1])

        return chain_tuple

    def get_macrocycle_marginal_atoms(self, file):
        try:
            macrofile = open(file, "r").readlines()

            return macrofile[0].split(" ")[0], macrofile[-1].split(" ")[0]
        except IOError:
            self.raise_popup_menu("An error occurred durning generation of file. Make sure you provided correct data.")

    def generate_bondfile(self):
        print "  Composing a bondfile from given component..."
        output_file = open(self._full_path_to_dir + os.sep + "bondfile.txt", "w")
        macro_names = []
        reg = re.compile("(\d)*_[\w]_[\d]*")

        for idx, component in enumerate(self.macrolink_list):
            if not all(len(loop[0].get()) == 0 and len(loop[1].get()) == 0 for loop in component):
                output_file.write("#" + str(idx + 1) + "\n")
                macro_name = ""
                for id, loop in enumerate(component):
                    beg = loop[0].get()
                    end = loop[1].get()
                    if not reg.match(beg) or not reg.match(end):
                        self.raise_popup_menu("Component " + str(idx + 1) + " line " + str(id + 1) +
                                              " contains invalid characters.")

                    if len(beg) == 0 and len(end) == 0:
                        break
                    if len(beg) < 5 or len(end) < 5:
                        output_file.close()
                        os.remove("bondfile.txt")
                        self.raise_popup_menu("There are missing data in the fields.")

                    first = beg.split("_")
                    second = end.split("_")
                    output_file.write(first[0] + "_" + first[1] + " " + first[2] + " " +
                                      second[0] + "_" + second[1] + " " + second[2] + "\n")
                    macro_name += (first[0] + first[1] + "_")
                macro_names.append(macro_name[:-1])
        output_file.close()

        return macro_names

    def convert_bondfile_to_components(self):
        print "  Generating component file..."
        macrolink = self._filenames[0]
        command = self.prog_converter + self._full_path_to_dir + os.sep + macrolink + " -x " + self._full_path_to_dir \
                  + os.sep + "bondfile.txt"
        res = os.popen(command).read()
        if "ERROR!!!" in res:
            err_msg = res.split("ERROR!!!")[1]
            self.raise_popup_menu(err_msg)

        if self.is_correct_data.get() and "WARNING!!!" in res:
            self.raise_popup_menu("Given structure has gaps or there are unnatural distances between atoms. "
                                  "Calculations aborted.")

    def rename_macrocomponent_files(self, mc):
        chain_tuple = []

        for idx, name in enumerate(mc):
            os.rename(system_working_directory + os.sep + "component_" + str(idx + 1) + ".pdb", "macro_" + name + ".pdb")
            os.rename(system_working_directory + os.sep + "component_" + str(idx + 1) + ".xyz", "macro_" + name + ".xyz")

            atoms = open(system_working_directory + os.sep + "macro_" + name + ".xyz").readlines()
            chain_tuple.append("macro_" + name + " " + atoms[0].split(" ")[0] + " " + atoms[-1].split(" ")[0])

        return chain_tuple

    ####################################################################################################################
    #                                               EXECUTE PROGRAM 2/2
    ####################################################################################################################

    def create_main_link_directory(self):
        self.current_working_dir = system_working_directory
        dir_name = self._filenames[0][:-3]
        directory = self.create_polymer_directory(dir_name[:-1])
        self.link_directory = self._full_path_to_dir + os.sep + directory

    def create_polymer_directory(self, prot):
        direct = self._full_path_to_dir + os.sep + prot

        if not os.path.exists(direct):
            os.makedirs(direct)

        return direct.split(os.sep)[-1]

    def move_files_to_polymer_directory(self, link_dir, hc):
        os.chdir(self.link_directory)
        self.separate_files_to_directory(self.current_working_dir, link_dir, hc + "_(inters(_OUT)*|results|glns).txt")
        self.separate_files_to_directory(self.current_working_dir, link_dir + os.sep + "_sm", hc + "_.*_sm\d*.xyz")
        self.separate_files_to_directory(self.current_working_dir, link_dir, hc + "_.*_cl\d*.jm")
        os.chdir(self.current_working_dir)

    def separate_files_to_directory(self, source, target, regex):
        if not os.path.exists(target):
            os.makedirs(target)

        reg = re.compile(regex)
        for filename in os.listdir(source):
            if reg.match(filename):
                shutil.move(os.path.join(source, filename), os.path.join(target, filename))

    def find_links(self):
        '''
            Since Windows does not provide a fork function crucial for multiprocessing, multithreading is not available.
            Therefore, each link needs to be checked in iteration which unfortunately prolongs the calculations.
            Thanks Bill Gates!
        '''
        print "  Calculating links..."
        self.output_data = {}
        self.selected_chains_links = {}
        self.all_complex_links = {}
        for idx, elem in enumerate(self.user_data):
            args = _calculate_links(idx, elem, self.output_data, self.selected_chains_links, self.all_complex_links)
            if args:
                self.output_data = self.merge_two_dicts(self.output_data, args[0])
                self.selected_chains_links = self.merge_two_dicts(self.selected_chains_links, args[1])
                self.all_complex_links = self.merge_two_dicts(self.all_complex_links, args[2])

    def merge_two_dicts(self, x, y):
        z = x.copy()
        z.update(y)
        return z

    def get_links_intersections(self):
        self.intersections = {}
        self.probabilistic_surface_idx = {}

        for comb in self.selected_chains_links.keys():
            comb = self.selected_chains_links[comb]
            if self.is_macrolink:
                chains = ["(" + str(_ + 1) + ")" for _ in range(len(comb))]
            elif all(".xyz" in loop for loop in comb):  # ions
                chains = [_.split(".pdb_")[1].split("_")[0] for _ in comb]
            else:
                chains = [_.split(" ")[0][-1] for _ in comb]
            hashcode = self.output_data[comb][1]

            for idx, link_name in enumerate(self.output_data[comb][0][1:]):
                link_name = link_name.split("x")[0]
                self.read_intersections(hashcode, link_name, idx, chains)

    def read_intersections(self, hc, l_n, arr_idx, ch):
        path_to_intersections = self.current_working_dir + os.sep + hc + "_inters_OUT.txt"

        with open(path_to_intersections, "r") as input:
            for idx, line in enumerate(input):
                line = line.rstrip().split("X")
                if line[0][:-2] == l_n:
                    retrieved_inters = line[1:]
                    pos = 0
                    intersections = ["" for _ in range(len(ch))]

                    for is_crossed in xrange(len(ch)):
                        for crossing in xrange(is_crossed + 1, len(ch)):
                            for cr in retrieved_inters[pos].split(" "):
                                if cr != "":
                                    intersections[is_crossed] += cr + ch[crossing] + ", "
                            for cr in retrieved_inters[pos + 1].split(" "):
                                if cr != "":
                                    intersections[crossing] += cr + ch[is_crossed] + ", "
                            pos += 2
                    if self.notebook.getcurselection() == "Probabilistic":
                        self.probabilistic_surface_idx[hc + str(arr_idx)] = str(idx)
                    self.intersections[hc + str(arr_idx)] = intersections
                    break

    def check_smoothing_avaliability(self):
        self.can_be_smoothed = {}

        if self.notebook.getcurselection() == "Probabilistic":
            self.calculate_probabilistic_smoothing()
        else:
            self.calculate_deterministic_smoothing()

    def calculate_probabilistic_smoothing(self):
        approach = " -2"
        num_tries = " 0"

        for cmb in self.selected_chains_links:
            comb = self.selected_chains_links[cmb]
            hashcode = self.output_data[comb][1]
            orig_links = self.output_data[comb][0][1:]
            for pos, elem in enumerate(orig_links):
                if int(float(elem.split("x")[1])) < 10:
                    break
            if pos == 0 or pos == len(orig_links)-1:
                orig_links = orig_links
            else:
                orig_links = orig_links[:pos]
            num_comp = str(len(comb))
            command = self.prog_homfly + num_comp + approach + num_tries + " "
            self.purge_lmknots(hashcode + "_")
            for idx, br in enumerate(comb):
                file = br.split(" ")[0][:-1]
                chain = br.split(" ")[0][-1]
                br_beg = br.split(" ")[1]
                br_end = br.split(" ")[2]
                comp = str(idx)
                path_to_dir = self.link_directory + os.sep + hashcode + os.sep + "_sm"
                smooth_filename = path_to_dir + os.sep + hashcode + "_ch" + comp + "_" + file + ".pdb" + \
                                  "_" + chain + "_L" + br_beg + "-" + br_end + "_clCLLVL_smSMLVL.xyz"
                smooth_ADD_filename = path_to_dir + os.sep + hashcode + "_ch" + comp + "ADD_" + file + ".pdb" + \
                                      "_" + chain + "_L" + br_beg + "-" + br_end + "_clCLLVL_smSMLVL.xyz"
                command += smooth_filename + " " + smooth_ADD_filename + " "
            command += "-ht " + hashcode

            for idx, topology in enumerate(orig_links):
                topology = topology.split("x")[0]
                closure = self.probabilistic_surface_idx[hashcode + str(idx)]
                smooth_lvl = self.smooth_level.get()

                while smooth_lvl > 0:
                    new_command = re.sub("SMLVL", str(smooth_lvl), command)
                    new_command = re.sub("CLLVL", str(closure), new_command)

                    call_homfly(new_command)
                    call_poly(hashcode)

                    try:
                        smooth_links = call_ncuc(hashcode, num_comp).split(" ")[3:-1]
                        smooth_links = self.fix_link_names(smooth_links)[1].split("x")[0].replace("X", " ")
                        self.purge_lmknots(hashcode + "_")
                    except IndexError:
                        print "crucial error with ncucLinks"
                        self.purge_main_directory()
                    if topology == smooth_links:
                        self.purge_lmknots(hashcode + "_")
                        self.can_be_smoothed[cmb + str(idx)] = [True, smooth_lvl]
                        break
                    else:
                        smooth_lvl = int(floor(float(smooth_lvl) / 2))
                else:
                    self.can_be_smoothed[cmb + str(idx)] = [False, None]
                    print " WARNING! Change of topology for ", file + chain, " topology", topology, \
                        ". Smoothing will not be available."

    def calculate_deterministic_smoothing(self):
        approach = " 0"
        num_tries = ""
        ions_in_links = False

        for cmb in self.output_data:
            hashcode = self.output_data[cmb][1]
            orig_links = self.output_data[cmb][0]
            num_comp = str(len(cmb))
            command = self.prog_homfly + num_comp + approach + num_tries + " 0 "
            self.purge_lmknots(hashcode + "_")
            filechains = []

            for idx, br in enumerate(cmb):
                path_to_dir = self.link_directory + os.sep + hashcode + os.sep + "_sm"
                comp = str(idx)
                closure = "0"
                br = br.split(" ")
                if self.is_macrolink:
                    file = br[0]
                    chain = ""
                    filechains.append(file)
                    smooth_filename = path_to_dir + os.sep + hashcode + "_ch" + comp + "_" + file + "_L"
                elif ".xyz" in br[0]:
                    ions_in_links = True
                    file = br[0][:-4]
                    chain = br[0].split(".pdb_")[1].split("_")[0]
                    filename = br[0].split(".pdb_")[0]
                    filechains.append(filename + chain)
                    smooth_filename = path_to_dir + os.sep + hashcode + "_ch" + comp + "_" + file + "_L"
                else:
                    file = br[0][:-1]
                    chain = br[0][-1]
                    filechains.append(file + chain)
                    smooth_filename = path_to_dir + os.sep + hashcode + "_ch" + comp + "_" + file + ".pdb" + "_" + \
                                      chain + "_L"
                br_beg = br[1]
                br_end = br[2]
                smooth_filename += br_beg + "-" + br_end + "_cl" + closure + "_smSMLVL.xyz"
                command += smooth_filename + " " + br_beg + " " + br_end + " "
            command += "-ht " + hashcode
            if ions_in_links:
                command += " -noCheckIds 1"

            smooth_lvl = self.smooth_level.get()
            while smooth_lvl > 0:
                new_command = re.sub("SMLVL", str(smooth_lvl), command)

                call_homfly(new_command)
                call_poly(hashcode)
                smooth_links = call_ncuc(hashcode, num_comp).split(" ")[3:-1]
                # oh Wanda, Wanda...
                if len(smooth_links) > 1 and "-" == smooth_links[0]:
                    smooth_links = smooth_links[-2:]
                smooth_links = self.fix_link_names(smooth_links)
                if orig_links == smooth_links:
                    self.purge_lmknots(hashcode + "_")
                    if self.is_macrolink or self.is_ions:
                        for elem in filechains:
                            self.can_be_smoothed[elem] = [True, smooth_lvl]
                    else:
                        self.can_be_smoothed[filechains[0]] = [True, smooth_lvl]
                    break
                smooth_lvl = int(floor(float(smooth_lvl) / 2))
                self.purge_lmknots(hashcode + "_")
            else:
                if self.is_macrolink or self.is_ions:
                    for elem in filechains:
                        self.can_be_smoothed[elem] = [False, None]
                else:
                    self.can_be_smoothed[filechains[0]] = [False, None]
                print " WARNING! Change of topology for ", file + chain, ". Smoothing will not be available"

    def purge_lmknots(self, hc=""):
        if hc + "EMCode.txt" in os.listdir(os.getcwd()):
            os.remove(hc + "EMCode.txt")
        if hc + "LMKNOT.txt" in os.listdir(os.getcwd()):
            os.remove(hc + "LMKNOT.txt")

    def fix_link_names(self, links):
        fixed_process = []
        for sele in links:
            elem = sele.split("(")
            l = elem[0] + "x" + elem[1][:-2]
            fixed_process.append(l)

        return fixed_process

    def calculate_smoothed_structures(self):
        if not all(val == False for val in self.can_be_smoothed.values()):
            ions_in_links = False

            for filechain in self.can_be_smoothed:
                if self.notebook.getcurselection() == "Probabilistic":
                    bridges = self.selected_chains_links[filechain[:-1]]
                    hashcode = self.output_data[self.selected_chains_links[filechain[:-1]]][1]
                    cl = self.probabilistic_surface_idx[hashcode + str(filechain[-1])]
                else:
                    bridges = self.selected_chains_links[filechain]
                    cl = "0"
                smooth_lvl = str(self.can_be_smoothed[filechain][1])
                num_comp = str(len(bridges))
                hashcode = self.output_data[bridges][1]
                path_to_smooth_dir = self.link_directory + os.sep + hashcode + os.sep + "_sm" + os.sep
                command = self.prog_homfly + num_comp + " 0 0 "
                for idx, bridge in enumerate(bridges):
                    bridges = bridge.split(" ")

                    if self.is_macrolink:
                        file = bridges[0]
                        smooth_file = hashcode + "_ch" + str(idx) + "_" + file + "_L"
                    elif ".xyz" in bridges[0]:
                        ions_in_links = True
                        file = bridges[0][:-4]
                        smooth_file = hashcode + "_ch" + str(idx) + "_" + file + "_L"
                    else:
                        file = bridges[0][:-1]
                        chain = bridges[0][-1]
                        smooth_file = hashcode + "_ch" + str(idx) + "_" + file + ".pdb_" + chain + "_L"
                    br_beg = bridges[1]
                    br_end = bridges[2]
                    smooth_file += br_beg + "-" + br_end + "_cl" + cl + "_sm" + smooth_lvl + ".xyz "
                    path_to_smooth_file = path_to_smooth_dir + smooth_file
                    command += path_to_smooth_file + br_beg + " " + br_end + " "
                if ions_in_links:
                    command += "-noCheckIds 1 "
                command += "-sf 2"
                call_homfly(command)
                self.purge_lmknots()
                self.separate_files_to_directory(self.current_working_dir, path_to_smooth_dir, "s\d*_" + hashcode +
                                                 "_ch\d*(ADD)*_")

    def separate_link_files_to_hashdirectories(self):
        print "  Separating hashfiles to proper directories..."

        for comb in self.output_data:
            curr_hash = self.output_data[comb][1]
            hash_link_dir = self.link_directory + os.sep + curr_hash
            self.move_files_to_polymer_directory(hash_link_dir, curr_hash)

    def purge_output_data(self):
        tmp_dict = {}
        for obj in self.selected_chains_links.keys():
            if self.selected_chains_links[obj] in self.output_data:
                tmp_dict[self.selected_chains_links[obj]] = self.output_data[self.selected_chains_links[obj]]
        self.output_data = tmp_dict

    def purge_main_directory(self):
        print "  Purging main directory..."
        reg = "(_extracted.pdb|_EMCode.txt|_LMKNOT.txt|_glns.txt|(_)*inters(_OUT)*.txt|_(s|ch)\d(ADD)*_(\d*\w*).pdb_\w_L(-)*\d*-(-)*\d*_cl\d*(_sm\d*)*.(jm|xyz))"

        for file in os.listdir(system_working_directory):
            regex = re.search(reg, file)
            if regex:
                try:
                    os.remove(file)
                except OSError:
                    pass

    ####################################################################################################################
    #                                          CREATE LINK WINDOW
    ####################################################################################################################

    def create_link_window(self):
        link_title = ""

        for f in self._filenames:
            link_title += f + ", "
        link_title = "[" + link_title[:-2] + "]"
        self.win_link_info = Pmw.Dialog(self.parent, title='Link information window ' + link_title + '',
                                        buttons=["Show links", "Exit"],
                                        command=self.invoke_link_window_buttons)
        self.win_link_info.resizable(0, 0)
        self.show_links()

    def show_links(self):
        self.win_link_info.invoke("Show links")
        self.win_link_info.component("buttonbox").delete(0)

    def invoke_link_window_buttons(self, clicked_btn):
        if clicked_btn == "Show links":
            self.window_parent = self.win_link_info.interior()
            self.display_protein_information_window()
        else:
            self.win_link_info.withdraw()
            self.win_link_info.destroy()

    def display_protein_information_window(self):
        self.link_inf_group_chains = Pmw.Group(self.window_parent, tag_text='Chains')
        self.link_inf_group_chains.grid(column=0, row=0, padx=5, pady=5)
        self.protein_array_btns = []

        if self.selected_chains_links:
            self.create_protein_row_names()
            self.create_protein_column()
            self.create_protein_rows()
            self.remove_trivial_directories()
            if not self.is_macrolink and self.all_complex_links and len(self.all_complex_links) > 1:
                self.create_all_links_frame()
            if self.is_macrolink:
                self.create_macroexplanation_interior()
            if len(self.protein_array_btns) == 1:
                self.protein_array_btns[0].invoke()
            self.create_protein_window_hints()
        else:
            self.protein_array_btns.append(tk.Text(self.link_inf_group_chains.interior(), bg="white", height=2 + 1,
                                                   padx=0, pady=0, wrap="word", bd=0, highlightthickness=0))
            self.configure_tkinter_text(self.protein_array_btns[-1], "There was a problem in parsing the request. "
                                                                     "Probably the loops selected have a common part, "
                                                                     "or there is an error in the structure uploaded. "
                                                                     "Please check your selection and the structure.", 1)
            self.protein_array_btns[-1].grid(column=1, row=1, columnspan=4, rowspan=2)

    def configure_tkinter_text(self, tk_txt_elem, text, height):
        tk_txt_elem.configure(state="normal")
        tk_txt_elem.insert("1.0", str(text))
        tk_txt_elem.tag_add("center", 1.0, "end")
        tk_txt_elem.tag_configure("center", justify='center')
        tk_txt_elem.tag_add("spacing1", 1.0, "end")
        tk_txt_elem.tag_config("spacing1", spacing1=5 + height)
        tk_txt_elem.configure(state="disabled")

    def create_protein_row_names(self):
        rows = []

        for idx, i in enumerate(["Link type", "Probability %", "Method"]):
            rows.append(tk.Label(self.link_inf_group_chains.interior(), text=i, width=self.protein_row_name_width,
                                 font=self.row_name_font, justify="center"))
            rows[-1].grid(column=1 + idx, row=0)

    def create_protein_column(self):
        protein_col_list = [self.selected_chains_links.keys()[0]] if self.is_macrolink else self.selected_chains_links.keys()

        for idx, filechain in enumerate(sorted(protein_col_list)):
            btn_name = "Details" if self.is_macrolink else filechain

            self.protein_array_btns.append(tk.Button(self.link_inf_group_chains.interior(), text=btn_name,
                                                     justify="center", width=self.protein_btn_width,
                                                     font=self.row_name_font,
                                                     command=lambda x=filechain, y=idx: self.display_chain_information_window(x, y)))
            self.protein_array_btns[-1].grid(column=0, row=idx + 1)

    def create_protein_window_hints(self):
        hint_view = Pmw.Balloon(self.window_parent, relmouse="both")

        for i in self.protein_array_btns:
            hint_view.bind(i, textwrap.fill("Displays the chain's links.", self.hint_width))

    def create_protein_rows(self):
        self.protein_array = []
        height = 1

        if self.selected_chains_links:
            protein_row_list = [self.selected_chains_links.keys()[0]] if self.is_macrolink else self.selected_chains_links
            for filechain in protein_row_list:
                row_element = []
                comb = self.selected_chains_links[filechain]
                link = self.output_data[comb][0][0].split("x")
                link_name = link[0].capitalize()
                link_probability = link[1][:-1]
                self.displayed_num_of_comp = str(len(self.selected_chains_links[filechain]))
                link_method = "Probabilistic" if self.notebook.getcurselection() == "Probabilistic" else "Deterministic"
                link_type = self.create_link_type_interior(link_name, self.link_inf_group_chains.interior(), height)
                row_element.append(link_type)

                for idx, elem in enumerate([link_probability, link_method]):
                    row_element.append(tk.Text(self.link_inf_group_chains.interior(), bg="white", height=height + 1,
                                               width=self.protein_row_width, padx=0, pady=1, wrap="word", bd=0,
                                               highlightthickness=0))
                    self.configure_tkinter_text(row_element[-1], elem, 1)
                self.protein_array.append(row_element)

            row = 1
            for filechain in self.protein_array:
                for idx, elem in enumerate(filechain):
                    elem.grid(column=idx+1, row=row)
                row += 1

    def remove_trivial_directories(self):
        print "  Removing useless link directories..."
        hash_dirs = []
        for key in self.selected_chains_links:
            hash_dirs.append(self.output_data[self.selected_chains_links[key]][1])

        for hashdir in os.listdir(self.link_directory):
            if hashdir not in hash_dirs and len(hashdir) == 4:
                shutil.rmtree(self.link_directory + os.sep + hashdir)

    def create_all_links_frame(self):
        self.link_inf_group_all_links = Pmw.Group(self.window_parent, tag_text='All complex links found')
        self.link_inf_group_all_links.grid(column=1, row=0, padx=5, pady=5)

        self.protein_all_links = Pmw.ScrolledFrame(self.link_inf_group_all_links.interior(), hull_width=300, usehullsize=1,
                                                   hull_height=90 + (len(self.selected_chains_links.keys())*10))
        self.protein_all_links.grid(sticky='swen', column=0, row=0, padx=5, pady=5)

        comb_title = tk.Label(self.protein_all_links.interior(), text="Combination", width=20,
                              font=self.row_name_font, justify="center")
        comb_title.grid(column=0, row=0)
        link_title = tk.Label(self.protein_all_links.interior(), text="Link", width=10,
                              font = self.row_name_font, justify="center")
        link_title.grid(column=1, row=0)

        for idx, link in enumerate(self.all_complex_links.keys()):
            links = self.all_complex_links[link].split("x")[0]
            link = ", ".join(link)
            height = max(textwrap.fill(link, 8).count("\n"), textwrap.fill(links, 8).count("\n"))

            content = tk.Text(self.protein_all_links.interior(), bg="white", width=25, padx=0, pady=0, wrap="word",
                              bd=0, highlightthickness=0, height=height)
            self.configure_tkinter_text(content, link, 0)
            content.grid(row=idx + 2, column=0)
            content = tk.Text(self.protein_all_links.interior(), bg="white", width=15, padx=0, pady=0, wrap="word",
                              bd=0, highlightthickness=0, height=height)
            self.configure_tkinter_text(content, links, 0)
            content.grid(row=idx + 2, column=1)

    def create_macroexplanation_interior(self):
        macroexplanation = ""
        for idx, filechain in enumerate(sorted(self.output_data.keys())[0]):
            macroexplanation += "Component " + str(idx + 1) + " = " + filechain + "\n"
        macroexplanation = macroexplanation[:-1]
        macrolink_label = tk.Label(self.link_inf_group_chains.interior(), text=macroexplanation)
        macrolink_label.grid(column=0, row=len(self.protein_array_btns) + 2, columnspan=4, pady=2, padx=5)

    def create_link_type_interior(self, link_type, interior, height, width=20):
        elem = tk.Text(interior, bg="white", height=height + 1, width=20, padx=0, pady=1, wrap="word", bd=0,
                       highlightthickness=0)
        try:
            if link_type.capitalize() == "Complex" or link_type.capitalize() == "Trivial":
                raise Exception
            linkname = " ".join([word.capitalize() for word in link_type.split(" ")]).replace("sup", "SUP")
            link_type, link_sup, link_sub, link_img_name = self.simplify_link_img(linkname)[1:]
            image = tk.PhotoImage(file=os.path.join(plugin_path + os.sep + "_links_img", link_img_name)).subsample(6, 6)
            tmp_label_lasso = tk.Label(self.win_link_info.interior())  # without label, pics will disappear
            tmp_label_lasso.image = image
            link_img = tk.Label(elem, image=image, bg="white", height=25, width=25, bd=0)
            link_img.grid(column=0, row=0, padx=5)

            link_type_name = tk.Text(elem, bg="white", height=height, padx=0, pady=0, wrap="word",
                                     width=width - 9, bd=0, highlightthickness=0)
            self.configure_linkname_text(link_type_name, link_type, link_sup, link_sub, height)
            link_type_name.grid(column=1, row=0)
        except Exception:
            link_type_name = tk.Text(elem, bg="white", height=height + 1, padx=0, pady=0, wrap="word", width=width,
                                     bd=0, highlightthickness=0)
            self.configure_tkinter_text(link_type_name, link_type, 1)
            link_type_name.grid(column=1, row=0)
        elem.configure(state="disabled")

        return elem

    def simplify_link_img(self, name, ncomp=0):
        if name == 'Other' or name == 'PUSTY': return '$Other$', 'Other', [], [], 'other.gif'
        calculated_ncomp = 0
        split_sum = name.split(' U ')
        for k in range(len(split_sum)):
            split_sum[k] = re.sub(' H ', ' # ', split_sum[k])
            split_sum[k] = split_sum[k].split('#')
        latex = '$'
        name_string = ''
        superscript = []
        subscript = []
        figure_name = ''
        for split in split_sum:
            for connected in split:
                latex_add, name_string_add, superscript_add, subscript_add, figure_name_add, k = self.gen_link_name_part(
                    connected, len(name_string))
                latex += latex_add + ' \# '
                name_string += name_string_add + ' # '
                superscript += superscript_add
                subscript += subscript_add
                figure_name += figure_name_add + 'h'
                calculated_ncomp += int(k)
            latex = latex[:-3] + '\cup '
            name_string = name_string[:-2] + 'u '
            calculated_ncomp = calculated_ncomp + 1 - len(split)
            figure_name = figure_name[:-1]
        try:
            ncomp = int(ncomp)
        except:
            ncomp = calculated_ncomp
        ncomp = max(int(ncomp), calculated_ncomp)
        for k in range(calculated_ncomp, ncomp):
            latex += '0_1 \cup '
            name_string += '01 u '
            subscript.append(len(name_string) - 4)
        if set(name_string[:-3].split()) == set(['Unlink', 'u', '01']):
            name_string = 'Unlink   '
            latex = '$Unlink      '
            subscript = []
            superscript = []
        return latex[:-6] + '$', name_string[:-3], superscript, subscript, figure_name + '_' + str(ncomp) + '.gif'

    def gen_link_name_part(self, part, total_len):
        part = part.strip()
        latex = ''
        name_string = ''
        superscript = []
        subscript = []
        figure_name = ''
        k = 0
        if part[0] == '+' or part[0] == '-' or part in ['3_1', '4_1', '6_3', '7_1', '8_3']:
            latex = part.split('_')[0] + '_{' + part.split('_')[1] + '}'
            name_string = re.sub('_', '', part)
            subscript.append(total_len + len(name_string) - 1)
            figure_name = re.sub('[-+]', '', name_string)
            k = 1
        elif part.strip() == '4-1':
            latex = '4_1'
            name_string = '41'
            subscript.append(total_len + 1)
            figure_name = '41'
            k = 1
        elif part[1] == '-':
            name, n = part.split('_')
            latex = name[0] + '_{' + name.split('-')[1] + '}'
            name_string = re.sub('-', '', name)
            figure_name = name_string
            for k in range(1, len(name_string)):
                subscript.append(total_len + k)
            for k in range(1, int(n)):
                latex += ' \cup 0_1'
                name_string += ' u 01'
                subscript.append(total_len + len(name_string) - 1)
            k = int(n)
        elif part[1:4] == 'SUP' and part.split('.')[0] != '2SUP2_1':
            latex = part[0] + '^' + part[4] + '_{' + part.split('_')[1].split('.')[0] + '}'
            if len(part.split('.')) > 1: latex += '.' + part.split('.')[1]
            name_string = re.sub('SUP', '', re.sub('_', '', part))
            figure_name = name_string.split('.')[0]
            superscript.append(total_len + 1)
            for k in range(2, len(name_string.split('.')[0])): subscript.append(total_len + k)
            k = int(part[4])
        elif part == '633_3.png':
            latex = '6^3_3'
            name_string = '633'
            superscript.append(total_len + 1)
            subscript.append(total_len + 2)
            figure_name = name_string
            k = 3
        elif 'Ring' in part:
            if part == 'Ring':
                latex = '0_1'
                k = 1
                name_string = '01'
                subscript.append(total_len + 1)
            else:
                k = int(part.split()[0])
                for n in range(k):
                    latex += '0_1 \cup '
                    name_string += '01 u '
                    subscript.append(total_len + len(name_string) - 4)
                latex = latex[:-6]
                name_string = name_string[:-3]
        elif 'Unlink' in part:
            k = int(part.split('_')[1])
            name_string = 'Unlink'
            latex = 'Unlink'
            figure_name = 'unlink'
        elif part[0] == 'L':
            latex = part
            name_string = part
            figure_name = part.split('.')[0]
            k = 4
        part = re.sub('_', '.', re.sub('2SUP2_1', 'Hopf', part))  # to deal with Hopf_2 and 2SUP2_1 ...
        if part.split('.')[0] in ['Hopf', 'Solomon', 'Whitehead', 'Star of David']:
            figure_name = re.sub('Star of David', 'davidstar', part.split('.')[0]).lower()
            k = 2
            part = re.sub('Star of David', '6^2_1', part)
            latex = part
            name_string = re.sub('[\^_]', '', part)
            if '^' in part:
                subscript.append(total_len + 2)
                superscript.append(total_len + 1)
        return latex, name_string, superscript, subscript, figure_name, k

    def configure_linkname_text(self, tk_txt_elem, content, sup, sub, h):
        tk_txt_elem.configure(state="normal")
        tk_txt_elem.tag_configure("subscript", offset=-5)
        tk_txt_elem.tag_configure("supscript", offset=5)

        for idx, elem in enumerate(content):
            if str(elem).isalpha() or str(elem).isdigit():
                if idx in sup:
                    tk_txt_elem.insert("insert", "", "", elem, "supscript")
                elif idx in sub:
                    tk_txt_elem.insert("insert", "", "", elem, "subscript")
                else:
                    tk_txt_elem.insert("insert", elem)
            else:
                tk_txt_elem.insert("insert", elem)

        tk_txt_elem.tag_add("center", 1.0, "end")
        tk_txt_elem.tag_configure("center", justify='center')
        tk_txt_elem.tag_add("spacing1", 1.0, "end")
        tk_txt_elem.tag_config("spacing1", spacing1=6 + h)
        tk_txt_elem.configure(state="disabled")

    def display_chain_information_window(self, filechain, idx):
        self.displayed_filechain = filechain
        self.displayed_link = idx
        self.color_line_in_array_of_results(self.protein_array, self.displayed_link)
        self.purge_pymol_structures()

        if hasattr(self, "link_inf_group_selected_chains_links") and self.link_inf_group_selected_chains_links.winfo_exists():
            self.link_inf_group_selected_chains_links.grid_forget()
        if hasattr(self, "gln_matrices") and self.gln_matrices.winfo_exists():
            self.gln_matrices.grid_forget()

        self.link_inf_group_selected_chains_links = Pmw.Group(self.window_parent,
                                                              tag_text='Links in file ' + str(filechain[:-1]) +
                                                                       ' and chain ' + str(filechain[-1]))
        self.link_inf_group_selected_chains_links.grid(column=0, row=1, padx=5, pady=5)

        self.create_chain_row_names()
        self.create_chain_rows()
        self.create_view_details_buttons()

        if self.is_gln_selected.get():
            self.display_gln_content()

        if not self.is_macrolink:
            self.pymol_display_chain(self.selected_chains_links[self.displayed_filechain])

        self.view_details_buttons[0].invoke()
        self.create_chart()
        self.create_display_links_hints()

    def create_chain_row_names(self):
        rows = []
        chains = self.selected_chains_links[self.displayed_filechain]
        row_names = ["Link type", "Probability %", "Ranges"]

        if all(chains[0].split(" ")[0][-1] == chain.split(" ")[0][-1] for chain in chains[1:]):
            if self.is_ions:
                for idx, loop in enumerate(chains):
                    row_names.append("Loop " + str(idx + 1))
                self.chain_array_row_name_width = int(79 / len(row_names))
                self.chain_array_row_width = self.chain_array_row_name_width + 8
            else:
                row_names += ["N loop piercings", "C loop piercings"]
                self.chain_array_row_name_width = 15
                self.chain_array_row_width = 23

            for idx, i in enumerate(row_names):
                rows.append(tk.Label(self.link_inf_group_selected_chains_links.interior(), text=i,
                                     width=self.chain_array_row_name_width,
                                     font=self.row_name_font, justify="center"))
                rows[-1].grid(column=1 + idx, row=0)
        else:
            if len([_.split(" ")[0][-1] for _ in chains]) == len(set([_.split(" ")[0][-1] for _ in chains])):
                for idx, ch in enumerate(chains):
                    if self.is_macrolink:
                        row_names.append("Component " + str(idx + 1))
                    else:
                        ch = ch.split(" ")[0][-1]
                        row_names.append("Chain " + ch + " piercings")
            else:
                for idx, ch in enumerate(chains):
                    if self.is_macrolink:
                        row_names.append("Component " + str(idx + 1))
                    else:
                        loop = ch.split(" ")[1] + "-" + ch.split(" ")[2]
                        row_names.append("Loop " + loop + " piercings")
            self.chain_array_row_name_width = int(84/len(row_names))
            self.chain_array_row_width = self.chain_array_row_name_width + 8

            for idx, i in enumerate(row_names):
                rows.append(tk.Label(self.link_inf_group_selected_chains_links.interior(), text=i,
                                     width=self.chain_array_row_name_width, font=self.row_name_font,
                                     justify="center"))
                rows[-1].grid(column=1 + idx, row=0)
            if len(chains) > 2:
                rows[0].configure(width=8)
                rows[1].configure(width=10)
                rows[2].configure(width=8)
                self.chain_array_row_width = [14, 14, 15]
                for elem in rows[3:]:
                    elem.configure(width=50/len(rows[3:]))
                    self.chain_array_row_width.append(68/len(rows[3:]))

    def color_line_in_array_of_results(self, array, selected):
        for idx, elem in enumerate(array):
            if idx == selected:
                for i in elem:
                    if len(i.winfo_children()):
                        for j in i.winfo_children():
                            j.configure(bg="lightblue")
                    i.configure(bg="lightblue")
            else:
                for i in elem:
                    if len(i.winfo_children()):
                        for j in i.winfo_children():
                            j.configure(bg="white")
                    i.configure(bg="white")

    def create_chain_rows(self):
        self.chain_array_rows = []
        sel_chain = self.selected_chains_links[self.displayed_filechain]

        for idx, link_elem in enumerate(self.output_data[sel_chain][0][1:]):
            self.displayed_num_of_comp = str(len(sel_chain))
            link_name = link_elem.split("x")[0]
            link_probability = link_elem.split("x")[1]
            if int(float(link_probability)) < 10:
                idx=-1
                break

            loop_ranges = ""
            for i, elem in enumerate(sel_chain):
                elem = elem.split(" ")

                if self.is_macrolink:
                    loop_ranges += elem[1] + "-" + elem[2] + "(" + str(i + 1) + ")" + ", "
                elif ".xyz" in elem[0]:  # ions
                    chain = elem[0].split(".pdb_")[1].split("_")[0]
                    loop_ranges += elem[1] + "-" + elem[2] + chain + ", "
                else:
                    loop_ranges += elem[1] + "-" + elem[2] + elem[0][-1] + ", "
            loop_ranges = loop_ranges[:-2]

            hashcode = str(self.output_data[sel_chain][1])
            intersections = self.intersections[hashcode + str(idx)]
            n_crossings = intersections[0]
            c_crossings = intersections[1]
            cell_height = max(textwrap.fill(loop_ranges, 16).count("\n"), textwrap.fill(n_crossings, 16).count("\n"),
                              textwrap.fill(c_crossings, 16).count("\n"), textwrap.fill(link_name, 16).count("\n")) + 2
            link_type = self.create_link_type_interior(link_name.lower(), self.link_inf_group_selected_chains_links.interior(),
                                                       cell_height, 25)

            row_values = [link_type, link_probability, loop_ranges]
            for elem in range(len(sel_chain)):
                try:
                    row_values.append(intersections[elem])
                except IndexError:
                    row_values.append("")

            row = []
            for num, elem in enumerate(row_values):
                width = self.chain_array_row_width if type(self.chain_array_row_width) == int else self.chain_array_row_width[num]
                if num != 0:
                    row.append(tk.Text(self.link_inf_group_selected_chains_links.interior(), bg="white", bd=0,
                                       height=cell_height, wrap="word", highlightthickness=0, width=width))
                    self.configure_tkinter_text(row[-1], elem, cell_height)
                else:
                    row.append(elem)
                row[-1].grid(column=num+1, row=idx + 1)
            self.chain_array_rows.append(row)

    def hasNumbers(self, inputString):
        return any(char.isdigit() for char in inputString)

    def create_view_details_buttons(self):
        self.view_details_buttons = []

        for idx in range(len(self.chain_array_rows)):
            self.view_details_buttons.append(tk.Button(self.link_inf_group_selected_chains_links.interior(),
                                                       text="view details", width=self.view_details_btn_width,
                                                       command=lambda x=idx: self.pymol_display_link(x)))
            self.view_details_buttons[-1].grid(sticky='w', column=0, row=idx+1)

    def create_chart(self):
        if (hasattr(self, "link_inf_probabilistic_chart") and self.link_inf_probabilistic_chart.winfo_exists()):
            self.link_inf_probabilistic_chart.grid_forget()

        self.link_inf_probabilistic_chart = Pmw.Group(self.window_parent, tag_text='Probability chart')
        self.link_inf_probabilistic_chart.grid(column=1, row=1, rowspan=len(self.protein_array_btns)+3, padx=5, pady=5)

        links = self.output_data[self.selected_chains_links[self.displayed_filechain]][0][1:]
        link_types = []
        link_probabilities = []
        colors = []
        others = 0

        for data in links:
            data = data.split("x")
            if int(float(data[1])) > 1:
                link_types.append(data[0])
                link_probabilities.append(float(data[1]))
            else:
                others += float(data[1])
        else:
            if sum(link_probabilities) < 100:
                if "Other" in link_types:
                    link_probabilities[link_types.index("Other")] = link_probabilities[link_types.index("Other")] + \
                                                                    (100 - sum(link_probabilities))
                else:
                    link_types.append("Other")
                    link_probabilities.append(100 - sum(link_probabilities))

        cmap = mplt.cm.get_cmap('hsv', len(link_types) + 1)
        for i in range(len(link_types)):
            colors.append(cmap(i))
        random.shuffle(colors)

        self.link_inf_chart_probability = mplt.figure.Figure(figsize=(self.chart_size, self.chart_size), facecolor='white')
        self.link_inf_chart_probability.subplots_adjust(top=0.99, bottom=0.01, left=0.05, right=0.97)
        subplot = self.link_inf_chart_probability.add_subplot(1, 1, 1)
        try:
            self.link_inf_piechart = subplot.pie(link_probabilities, labels=link_types, autopct='%.1f% %', colors=colors,
                                                 shadow=False, wedgeprops={"linewidth": 0})
        except TypeError:
            self.link_inf_piechart = subplot.pie(link_probabilities, labels=link_types, autopct='%.1f% %',
                                                 colors=colors, shadow=False)
        subplot.axis('equal')
        [i.set_visible(False) for i in self.link_inf_piechart[1]]
        [i.set_visible(False) for i in self.link_inf_piechart[2]]

        centre_circle = pyplt.Circle((0, 0), 0.67, color='black', fc='white', linewidth=0)
        ax = self.link_inf_chart_probability.add_subplot(1, 1, 1)
        self.link_inf_piechart_subtitle = ax.text(0, -0.55, "", ha="center", size=self.chart_subtitle)
        mplt.rcParams.update({'font.size': self.chart_font_size})
        ax.add_patch(centre_circle)
        mplt.rcParams.update(mplt.rcParamsDefault)

        self.link_inf_probability_canvas = FigureCanvasTkAgg(self.link_inf_chart_probability,
                                                             master=self.link_inf_probabilistic_chart.interior())
        self.link_inf_probability_canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        self.link_inf_probability_canvas.show()
        self.link_inf_probability_canvas.mpl_connect('motion_notify_event', self.display_chosen_probability)

        self.link_inf_chart_probability.add_axes(self.chart_img_axes)
        self.link_inf_chart_probability.axes[1].get_yaxis().set_visible(False)
        self.link_inf_chart_probability.axes[1].get_xaxis().set_visible(False)
        self.link_inf_chart_probability.axes[1].spines["top"].set_visible(False)
        self.link_inf_chart_probability.axes[1].spines["bottom"].set_visible(False)
        self.link_inf_chart_probability.axes[1].spines["left"].set_visible(False)
        self.link_inf_chart_probability.axes[1].spines["right"].set_visible(False)
        self.link_inf_chart_probability.axes[1].patch.set_visible(False)

        idx_most_probable_link = link_probabilities.index(max(link_probabilities))
        self.initialise_probability_chart_image(idx_most_probable_link, link_types[idx_most_probable_link])

    def initialise_probability_chart_image(self, idx, link_type):
        link_probability = self.link_inf_piechart[2][idx].get_text()
        content = self.simplify_link_img(link_type)
        txt = content[0][0] + "\\rm " + content[0][1:]
        filename = plugin_path + os.sep + "_links_img" + os.sep + content[-1]
        self.link_inf_chart_probability.axes[1].images = []
        try:
            img = mpimg.imread(filename)
            self.link_inf_chart_probability.axes[1].imshow(img, zorder=0)
            self.link_inf_piechart_subtitle.set_text(txt + "\n" + link_probability)
        except IOError:
            self.link_inf_chart_probability.axes[1].images = []
        self.link_inf_probability_canvas.draw()

    def display_chosen_probability(self, event):
        for idx, i in enumerate(self.link_inf_piechart[0]):
            (hit, _) = i.contains(event)
            if hit:
                link_type = i.get_label()
                link_probability = self.link_inf_piechart[2][idx].get_text()
                content = self.simplify_link_img(link_type)
                txt = content[0][0] + "\\rm " + content[0][1:]
                filename = plugin_path + os.sep + "_links_img" + os.sep + content[-1]
                self.link_inf_chart_probability.axes[1].images = []
                try:
                    img = mpimg.imread(filename)
                    self.link_inf_chart_probability.axes[1].imshow(img, zorder=0)
                    self.link_inf_piechart_subtitle.set_text(txt + "\n" + link_probability)
                except IOError:
                    self.link_inf_chart_probability.axes[1].images = []
                self.link_inf_probability_canvas.draw()

    def create_display_links_hints(self):
        hint_view = Pmw.Balloon(self.window_parent, relmouse="both")

        for elem in self.view_details_buttons:
            hint_view.bind(elem, textwrap.fill("Displays the chosen link.", self.hint_width))

    def pymol_display_link(self, idx):
        self.displayed_link = idx
        comb = self.selected_chains_links[self.displayed_filechain]
        self.color_line_in_array_of_results(self.chain_array_rows, self.displayed_link)

        is_ion = False
        if not self.is_macrolink:
            self.pymol_display_chain(self.selected_chains_links[self.displayed_filechain])
        hashcode = self.output_data[comb][1]
        closure = str(self.probabilistic_surface_idx[hashcode + str(idx)]) if self.notebook.getcurselection() == "Probabilistic" else "0"

        for elem in cmd.get_names(type="all"):
            if elem.startswith("SURF") or elem.startswith("PIERC") or elem.startswith("BR") or elem.startswith("SM") \
                    or elem.startswith("BR") or elem.startswith("SM") or elem.startswith("macro_"):
                cmd.delete(name=elem)
        cmd.hide(representation="spheres", selection="all")
        cmd.hide(representation="nonbonded", selection="all") # hide nonbonded dots for macrolinks
        if hasattr(self, "ions_connections") and self.ions_connections:
            tmp_pos_list = [self.is_ion_link_clicked.index(i) for i in self.is_ion_link_clicked if i.get() == 1]

        if self.is_macrolink:
            for id, link in enumerate(comb):
                filename, file_beg, file_end = link.split(" ")
                path_to_surface = self.link_directory + os.sep + str(hashcode) + os.sep + hashcode + "_s" + str(id) \
                                  + "_" + filename + "_L" + file_beg + "-" + file_end + "_cl" + closure + ".jm"
                path_to_macrofile = self.link_directory + os.sep + filename + ".pdb"
                pymol_filename = filename

                cmd.load(filename=path_to_macrofile, object=pymol_filename)
                self.color_line_in_array_of_results(self.chain_array_rows, self.displayed_link)
                self.get_cgo_coordinates(path_to_surface)
                self.pymol_draw_piercings(id, pymol_filename)
                self.pymol_draw_surface(id, pymol_filename)
                cmd.deselect()
                cmd.center(selection=pymol_filename)
                cmd.orient(selection=pymol_filename)
                cmd.set(name="line_width", value="4", selection=pymol_filename)
                cmd.show(representation="sticks", selection=pymol_filename)
                cmd.color(color=self.colors_of_sequence[id], selection=pymol_filename)
        else:
            for id, link in enumerate(comb):
                link = link.split(" ")
                path_to_dir = self.link_directory + os.sep + str(hashcode)
                if ".xyz"in link[0]: # ions
                    is_ion = True
                    filename = link[0][:-4]
                    chain = filename.split(".pdb_")[1].split("_")[0]
                    pymol_filename = filename.split(".pdb_")[0]
                    position = tmp_pos_list[id]
                    path_to_surface = path_to_dir + os.sep + hashcode + "_s" + str(id) + "_" + filename + "_L"
                    atoms = self.ions_connections[position][2].split(" ")
                    self.pymol_display_ion_structures(atoms, False)
                else:
                    filename = link[0][:-1]
                    chain = link[0][-1]
                    pymol_filename = filename
                    path_to_surface = path_to_dir + os.sep + hashcode + "_s" + str(id) + "_" + filename + ".pdb_" + \
                                      chain + "_L"

                bridge = [link[2], link[1]] if int(link[1]) > int(link[2]) else [link[1], link[2]]
                file_beg = link[1]
                file_end = link[2]
                path_to_surface += file_beg + "-" + file_end + "_cl" + closure + ".jm"
                atom = "ca"
                bridges = bridge[0] + "_" + bridge[1]
                self.color_line_in_array_of_results(self.chain_array_rows, self.displayed_link)
                self.get_cgo_coordinates(path_to_surface)
                if is_ion:
                    seq = self.generate_ion_selection(self.ions_connections[position][2].split(" "))
                    cmd.select("SEQ_" + chain + "_" + bridge[0] + "_" + bridge[1],
                               selection=pymol_filename + " and " + seq)
                else:
                    cmd.select("SEQ_" + chain + "_" + bridge[0] + "_" + bridge[1],
                               selection=pymol_filename + " and chain " + chain + " and residue " + bridge[0]
                                         + "-" + bridge[1])
                self.pymol_draw_piercings(id, chain + "_" + bridges)
                self.pymol_draw_surface(id, chain + "_" + bridges)
                if self.notebook.getcurselection() == "Probabilistic" and self.br_cylinders:
                    self.pymol_draw_cylinders("BR_" + chain + "_" + bridges)
                elif not is_ion:
                    cmd.select("BR_" + chain + "_" + bridge[0] + "_" + bridge[1],
                               selection="(" + pymol_filename + " and chain " + chain + " and residue " + bridge[0] +
                                         " and name " + atom + ")+" + "(" + pymol_filename + " and chain " + chain +
                                         " and residue " + bridge[1] + " and name " + atom + ")")
                    cmd.bond(atom1=pymol_filename + " and chain " + chain + " and residue " + bridge[0] + " and name " + atom,
                             atom2=pymol_filename + " and chain " + chain + " and residue " + bridge[1] + " and name " + atom)
                    cmd.show(representation='sphere', selection="BR_" + chain + "_" + bridge[0] + "_" + bridge[1])
                    cmd.show(representation="sticks", selection="BR_" + chain + "_" + bridge[0] + "_" + bridge[1])
                cmd.color(color=self.colors_of_sequence[id], selection="SEQ_" + chain + "_" + bridge[0] + "_" + bridge[1])
                cmd.deselect()
            print "  Surfaces drawn in PyMOL..."
            cmd.center(selection=pymol_filename + "_" + chain)
            cmd.orient(selection=pymol_filename + "_" + chain)
        cmd.clip(mode="slab", distance="1000")
        cmd.set(name="two_sided_lighting", value=1)
        self.mark_crossings_on_structure()
        self.create_bottom_buttons()

    def get_cgo_coordinates(self, path_to_file):
        if not os.path.isfile(path_to_file):
            self.raise_popup_menu('File with coordinates of vertices not found.')

        input_file = open(path_to_file, "r")
        self.piercing_coord = []
        self.surface_coord = []
        self.shallow_coord = []
        self.br_cylinders = []
        surface_coord = []
        polygon_num = ""

        for idx, line in enumerate(input_file.readlines()):
            fixed_line = list(re.findall(r"(draw (polygon|cylinder)\d(_intreduced)*)|(\\d|-?\d{1,3}[,.]\d{1,2})", line))
            if len(fixed_line) == 10:
                coord = list((re.findall(r"-?\d{1,3}[,.]\d{1,2}", line)))
                if "_intreduced" in fixed_line[0]:
                    self.shallow_coord.append(coord)
                elif polygon_num == "":
                    polygon_num = fixed_line[0]
                    surface_coord.append(coord)
                elif fixed_line[0] != polygon_num:
                    polygon_num = fixed_line[0]
                    crossing = surface_coord.pop()
                    self.piercing_coord = crossing
                    self.surface_coord.append(surface_coord)
                else:
                    surface_coord.append(coord)
            elif len(fixed_line) == 8:
                coord = list((re.findall(r"-?\d{1,3}[,.]\d{1,2}", line)))
                self.br_cylinders.append(coord)
        else:
            crossing = surface_coord.pop()
            self.piercing_coord = crossing
            self.surface_coord = surface_coord

    def generate_ion_selection(self, ion_path):
        path = ""
        for at1, rel, at2 in zip(ion_path[::2], ion_path[1::2], ion_path[2::2]):
            at1 = at1.split("_")
            at2 = at2.split("_")

            if int(at1[1][1:]) > int(at2[1][1:]): # if the first residue is bigger than second, we have to swap them
                tmp = at1
                at1 = at2
                at2 = tmp

            if rel == "...":
                path += "(chain " + at1[1][0] + " and residue " + at1[1][1:] + "-" + at2[1][1:] + " and name ca)+"

        return path[:-1]

    def mark_crossings_on_structure(self):
        piercings = []
        atom = "ca"
        pos_pierc = ""
        neg_pierc = ""
        for elem in self.chain_array_rows[self.displayed_link][3:]:
            piercings += str(elem.get("1.0", "end-1c")) + ", "
        piercings = filter(len, "".join(piercings).split(", "))

        if len(piercings) != 1:
            for pierc in piercings:
                if self.is_macrolink:
                    comp_num = int(pierc[-2]) - 1
                    chain = "X"
                    pierc = pierc[0] + str(int(pierc[1:-3]) * 2)
                    seq = 2
                    file = self.selected_chains_links[self.displayed_filechain][comp_num].split(" ")[0]
                else:
                    chain = pierc[-1]
                    pierc = pierc[:-1]
                    seq = 1
                    file = self.displayed_filechain[:-1]

                if "," in pierc:
                    if pierc.startswith("+"):
                        pos_pierc += "(" + file + " and chain " + chain + " and residue " + str(pierc[1:-1]) + \
                                     " and name " + atom + ") "
                        pos_pierc += "(" + file + " and chain " + chain + " and residue " + str(int(pierc[1:-1]) + seq)\
                                     + " and name " + atom + ") "
                    elif pierc.startswith("-"):
                        neg_pierc += "(" + file + " and chain " + chain + " and residue " + str(pierc[1:-1]) + \
                                     " and name " + atom + ") "
                        neg_pierc += "(" + file + " and chain " + chain + " and residue " + str(int(pierc[1:-1]) + seq)\
                                     + " and name " + atom + ") "
                else:
                    if pierc.startswith("+"):
                        pos_pierc += "(" + file + " and chain " + chain + " and residue " + str(pierc[1:]) + \
                                     " and name " + atom + ") "
                        pos_pierc += "(" + file + " and chain " + chain + " and residue " + str(int(pierc[1:]) + seq) \
                                     + " and name " + atom + ") "
                    elif pierc.startswith("-"):
                        neg_pierc += "(" + file + " and chain " + chain + " and residue " + str(pierc[1:]) + \
                                     " and name " + atom + ") "
                        neg_pierc += "(" + file + " and chain " + chain + " and residue " + str(int(pierc[1:]) + seq) \
                                     + " and name " + atom + ") "
        else:
            piercings = piercings[0]
            if piercings.endswith(","):
                piercings = piercings[:-1]
            if piercings.startswith("+"):
                chain = piercings[-1] if not self.is_macrolink else "X"
                piercings = piercings[:-1]
                pos_pierc += "(" + file + " and chain " + chain + " and residue " + str(piercings[1:]) + " and name " \
                             + atom + ") "
                pos_pierc += "(" + file + " and chain " + chain + " and residue " + str(
                    int(piercings[1:]) + 1) + " and name " + atom + ") "
            elif piercings.startswith("-"):
                chain = piercings[-1] if not self.is_macrolink else "X"
                piercings = piercings[:-1]
                neg_pierc += "(" + file + " and chain " + chain + " and residue " + str(piercings[1:]) + " and name " \
                             + atom + ") "
                neg_pierc += "(" + file + " and chain " + chain + " and residue " + str(
                    int(piercings[1:]) + 1) + " and name " + atom + ") "

        if len(pos_pierc) > 0:
            cmd.select(name="POS_PIERC", selection=pos_pierc[:-1])
            cmd.color(color="blue", selection="POS_PIERC")
        if len(neg_pierc) > 0:
            cmd.select(name="NEG_PIERC", selection=neg_pierc[:-1])
            cmd.color(color="green", selection="NEG_PIERC")
        cmd.deselect()

    def create_bottom_buttons(self):
        self.link_inf_surface_bridges = []

        self.link_inf_is_surface_displayed = [tk.IntVar()]
        self.link_inf_surface_button = tk.Checkbutton(self.link_inf_group_selected_chains_links.interior(),
                                                      text="Surface", anchor="center", indicatoron=0, width=10, height=2,
                                                      relief="ridge", variable=self.link_inf_is_surface_displayed[0],
                                                      command=lambda x=None, y=0: self.pymol_display_surface(x, y))
        self.link_inf_surface_button.select()
        self.link_inf_surface_button.grid(column=1, row=len(self.chain_array_rows) + 3, padx=10, pady=10)

        self.link_inf_is_smoothed = tk.IntVar()
        self.link_inf_smooth_button = tk.Checkbutton(self.link_inf_group_selected_chains_links.interior(),
                                                     text="Smooth", anchor="center", indicatoron=0, width=10, height=2,
                                                     relief="ridge", variable=self.link_inf_is_smoothed,
                                                     command=lambda x=self.displayed_filechain:
                                                     self.pymol_display_smooth(x))
        self.link_inf_smooth_button.grid(column=3, row=len(self.chain_array_rows) + 3, padx=10, pady=10)

        self.create_surfaces_button()
        self.create_surface_hints()

    def create_surfaces_button(self):
        key = self.selected_chains_links[self.displayed_filechain]
        hint_surfaces_button = Pmw.Balloon(self.window_parent, relmouse="both")

        for i in range(int(self.displayed_num_of_comp)):
            if self.is_macrolink:
                bridge_label, beg, end = key[i].split(" ")
                chain = "X"
                surf_text = "Surface (component " + str(i + 1) + ")"
            else:
                if int(key[i].split(" ")[1]) > int(key[i].split(" ")[2]):
                    beg = key[i].split(" ")[2]
                    end = key[i].split(" ")[1]
                else:
                    beg = key[i].split(" ")[1]
                    end = key[i].split(" ")[2]
                if ".xyz" in key[i].split(" ")[0]:
                    chain = key[i].split(" ")[0].split(".pdb_")[1].split("_")[0]
                    surf_text = "Surface (loop " + str(i + 1) + ")"
                else:
                    chain = key[i].split(" ")[0][-1]
                    surf_text = "Surface (" + beg + "-" + end + chain + ")"

                bridge_label = beg + "-" + end + chain
            surface_btn_width = self.surf_btn_width if len(surf_text) < self.surf_btn_width else len(surf_text) - 2

            self.link_inf_is_surface_displayed.append(tk.IntVar())
            self.link_inf_surface_bridges.append(tk.Checkbutton(self.link_inf_group_selected_chains_links.interior(),
                                                                text=surf_text, anchor="center", indicatoron=0,
                                                                relief="ridge", height=2, width=surface_btn_width,
                                                                variable=self.link_inf_is_surface_displayed[-1],
                                                                command=lambda x=bridge_label, y=i+1:
                                                                self.pymol_display_surface(x, y)))
            self.link_inf_surface_bridges[-1].grid(column=2*i, row=len(self.chain_array_rows) + 4, padx=10, columnspan=2,
                                                   pady=10)
            self.link_inf_surface_bridges[-1].select()
            hint_surfaces_button.bind(self.link_inf_surface_bridges[-1], textwrap.fill("Shows/Hides the surface for "
                                                                                       "bridge spanned on atoms " + beg
                                                                                       + " and " + end + " on chain " +
                                                                                       chain + ".", self.hint_width))

    def create_surface_hints(self):
        hint_bottom_btns = Pmw.Balloon(self.window_parent, relmouse="both")

        hint_bottom_btns.bind(self.link_inf_surface_button, "Shows/hides all surfaces spanned on the closed loop.")
        hint_bottom_btns.bind(self.link_inf_smooth_button, "Shows/hides the smoothed representation of the selected chain.")

    ####################################################################################################################
    #                                    METHODS OPERATING ON OBJECTS IN PYMOL
    ####################################################################################################################

    def pymol_draw_piercings(self, idx, pymol_cgo_name):
        color = self.colors_of_piercings[idx]
        triang = [BEGIN, TRIANGLES]

        n = self.compute_normal(float(self.piercing_coord[0]), float(self.piercing_coord[1]), float(self.piercing_coord[2]),
                                float(self.piercing_coord[3]), float(self.piercing_coord[4]), float(self.piercing_coord[5]),
                                float(self.piercing_coord[6]), float(self.piercing_coord[7]), float(self.piercing_coord[8]))
        triang.extend([
            COLOR, color[0], color[1], color[2],
            NORMAL, n[0], n[1], n[2],
            VERTEX, float(self.piercing_coord[0]), float(self.piercing_coord[1]), float(self.piercing_coord[2]),
            VERTEX, float(self.piercing_coord[3]), float(self.piercing_coord[4]), float(self.piercing_coord[5]),
            VERTEX, float(self.piercing_coord[6]), float(self.piercing_coord[7]), float(self.piercing_coord[8]),
        ])
        triang.append(END)
        cmd.load_cgo(triang, "PIERC_" + pymol_cgo_name)

    def compute_normal(self, x1, y1, z1, x2, y2, z2, x3, y3, z3):
        """
            Using coordinates this method computes normal perpendicular to the surface. Used to change the direction
            of light in PyMOL thus the surface of triangles is more visible.
        :return: vector containing normal perpendicular to the surface.
        """
        nx = (y2 - y1) * (z3 - z2) - (z2 - z1) * (y3 - y2)
        ny = (z2 - z1) * (x3 - x2) - (x2 - x1) * (z3 - z2)
        nz = (x2 - x1) * (y3 - y2) - (y2 - y1) * (x3 - x2)

        return nx, ny, nz

    def pymol_draw_surface(self, idx, pymol_cgo_name):
        color = self.colors_of_surfaces[idx]
        triang = [BEGIN, TRIANGLES]

        for vertex in self.shallow_coord:
            n = self.compute_normal(float(vertex[0]), float(vertex[1]), float(vertex[2]),
                                    float(vertex[3]), float(vertex[4]), float(vertex[5]),
                                    float(vertex[6]), float(vertex[7]), float(vertex[8]))
            triang.extend([
                COLOR, 0.8, 0.8, 0.8,
                NORMAL, n[0], n[1], n[2],
                VERTEX, float(vertex[0]), float(vertex[1]), float(vertex[2]),
                VERTEX, float(vertex[3]), float(vertex[4]), float(vertex[5]),
                VERTEX, float(vertex[6]), float(vertex[7]), float(vertex[8]),
            ])
        for vertex in self.surface_coord:
            n = self.compute_normal(float(vertex[0]), float(vertex[1]), float(vertex[2]),
                                    float(vertex[3]), float(vertex[4]), float(vertex[5]),
                                    float(vertex[6]), float(vertex[7]), float(vertex[8]))
            triang.extend([
                COLOR, color[0], color[1], color[2],
                NORMAL, n[0], n[1], n[2],
                VERTEX, float(vertex[0]), float(vertex[1]), float(vertex[2]),
                VERTEX, float(vertex[3]), float(vertex[4]), float(vertex[5]),
                VERTEX, float(vertex[6]), float(vertex[7]), float(vertex[8]),
            ])
        triang.append(END)
        cmd.load_cgo(triang, "SURF_" + pymol_cgo_name)

    def pymol_draw_cylinders(self, pymol_cgo_name):
        cylinder = []

        for vertex in self.br_cylinders:
            cylinder.extend([
                CYLINDER, float(vertex[0]), float(vertex[1]), float(vertex[2]),
                float(vertex[3]), float(vertex[4]), float(vertex[5]),
                0.5, 1, 0.5, 0, 1, 0.5, 0,
            ])
        cmd.load_cgo(cylinder, pymol_cgo_name)

    def pymol_display_surface(self, bridges, idx):
        if self.displayed_link is None:
            self.lasinf_surface_button.select()
            self.raise_popup_menu('No link chosen. Please load an appropriate link into the PyMOL viewer '
                                  '(e.g. press the button ''view details'') and try again.')

        if self.link_inf_is_smoothed.get():
            is_smoothed = "SM_"
            seq = "SEQ_SM_"
        else:
            is_smoothed = ""
            seq = "SEQ_"
        if not bridges:
            bridges = []
            id = 0
            for elem in self.selected_chains_links[self.displayed_filechain]:
                elem = elem.split(" ")
                if self.is_macrolink and "macro_" in elem[0]:  # macro
                    bridges.append(elem[0])
                elif ".xyz" in elem[0]:  # ions
                    br = [elem[2], elem[1]] if int(elem[1]) > int(elem[2]) else [elem[1], elem[2]]
                    if self.link_inf_is_smoothed.get():
                        chain = elem[0].split(" ")[0].split(".pdb_")[1].split("_")[0]
                        bridges.append(str(id) + "_" + br[0] + "-" + br[1] + chain)
                        id += 1
                    else:
                        chain = elem[0].split(" ")[0].split(".pdb_")[1].split("_")[0]
                        bridges.append(br[0] + "-" + br[1] + chain)
                else:
                    bridges.append(elem[1] + "-" + elem[2] + elem[0][-1])
        else:
            bridges = [bridges]

        if self.link_inf_is_surface_displayed[idx].get():
            for i, bridge in enumerate(list(bridges)):
                if self.is_macrolink:
                    cmd.show(representation="cgo", selection="SURF_" + is_smoothed + bridge)
                    cmd.show(representation="cgo", selection="PIERC_" + is_smoothed + bridge)
                else:
                    bridge = self.fix_bridges(bridge)
                    if self.is_ions and self.link_inf_is_smoothed.get():
                        if idx:
                            chain = str(idx - 1)
                            br_beg = bridge[0]
                            br_end = bridge[1][:-1]
                        else:
                            chain = bridge[0].split("_")[0]
                            br_beg = bridge[0].split("_")[1]
                            br_end = bridge[1][:-1]
                    else:
                        chain = bridge[1][-1]
                        br_beg = bridge[0]
                        br_end = bridge[1][:-1]

                    cmd.show(representation="cgo", selection="SURF_" + is_smoothed + chain + "_" + br_beg + "_" + br_end)
                    cmd.show(representation="cgo", selection="PIERC_" + is_smoothed + chain + "_" + br_beg + "_" + br_end)
                if idx == 0:
                    self.link_inf_surface_bridges[i].select()
                    if self.is_macrolink:
                        cmd.color(color=self.colors_of_sequence[i], selection=bridge)
                    else:
                        cmd.color(color=self.colors_of_sequence[i], selection=seq + chain + "_" + br_beg + "_" + br_end)
                else:
                    if self.is_macrolink:
                        cmd.color(color=self.colors_of_sequence[idx - 1], selection=bridge)
                    else:
                        cmd.color(color=self.colors_of_sequence[idx - 1], selection=seq + chain + "_" + br_beg + "_" + br_end)
                if all(elem.get() == 1 for elem in self.link_inf_is_surface_displayed[1:]):
                    self.link_inf_surface_button.select()
            self.mark_crossings_on_structure()
        else:
            for i, bridge in enumerate(list(bridges)):
                if self.is_macrolink:
                    bridge = bridge.split(" ")[0]

                    cmd.hide(representation="cgo", selection="SURF_" + is_smoothed + bridge)
                    cmd.hide(representation="cgo", selection="PIERC_" + is_smoothed + bridge)
                else:
                    bridge = self.fix_bridges(bridge)
                    if self.is_ions and self.link_inf_is_smoothed.get():
                        if idx:
                            chain = str(idx - 1)
                            br_beg = bridge[0]
                            br_end = bridge[1][:-1]
                            orig_chain = bridge[1][-1]
                            color = self.get_seq_color(orig_chain, br_beg, br_end)
                        else:
                            chain = bridge[0].split("_")[0]
                            br_beg = bridge[0].split("_")[1]
                            br_end = bridge[1][:-1]
                            orig_chain = bridge[1][-1]
                            color = self.get_seq_color(orig_chain, br_beg, br_end)
                    else:
                        chain = bridge[1][-1]
                        br_beg = bridge[0]
                        br_end = bridge[1][:-1]
                        color = self.get_seq_color(chain, br_beg, br_end)

                    cmd.hide(representation="cgo", selection="SURF_" + is_smoothed + chain + "_" + br_beg + "_" + br_end)
                    cmd.hide(representation="cgo", selection="PIERC_" + is_smoothed + chain + "_" + br_beg + "_" + br_end)
                    util.cbc(selection=seq + chain + "_" + br_beg + "_" + br_end, first_color=color, legacy=1)
                if idx == 0:
                    self.link_inf_surface_bridges[i].deselect()
            if any(elem.get() == 0 for elem in self.link_inf_is_surface_displayed[1:]):
                self.link_inf_surface_button.deselect()
            print "  Surface without area of triangulation showed..."

    def fix_bridges(self, br):
        """ Method to fix bridges in case there are negative indices"""
        bridge = br.split("-")
        if len(bridge) == 3:
            if bridge[0] == "":
                bridge = ["-" + bridge[1], bridge[2]]
            else:
                bridge = [bridge[0], "-" + bridge[1]]
        elif len(bridge) == 4:
            bridge = ["-" + bridge[1], "-" + bridge[3]]
        return bridge

    def get_seq_color(self, c, beg, end):
        for idx, elem in enumerate(self.selected_chains_links[self.displayed_filechain]):
            if ".xyz" in elem:  # ions
                elem = elem.split(" ")
                chain = elem[0].split(".pdb_")[1].split("_")[0]
                br = [elem[2], elem[1]] if int(elem[1]) > int(elem[2]) else [elem[1], elem[2]]
                filename = elem[0].split(" ")[0].split(".pdb_")[0] + ".pdb"
            else:
                if c + " " + beg + " " + end in elem:
                    filename = elem.split(" ")[0][:-1] + ".pdb"
        color = util._color_cycle[self._chains[self._filenames.index(filename)].index(str(c))]

        return color

    def pymol_display_smooth(self, filechain):
        if self.displayed_link is None:
            self.link_inf_smooth_button.deselect()
            self.raise_popup_menu('No link chosen. Please load an appropriate link into the PyMOL viewer '
                                  '(e.g. press the button ''view details'') and try again.')

        if (self.notebook.getcurselection() == "Probabilistic" and not self.can_be_smoothed[filechain +
            str(self.displayed_link)][0]) or (not self.notebook.getcurselection() == "Probabilistic" and not
        self.can_be_smoothed[filechain][0]):
            self.link_inf_smooth_button.deselect()
            self.raise_popup_menu("No smoothing is available due to a change in the topology.")

        if self.link_inf_is_surface_displayed[0].get() == 0:
            self.link_inf_surface_button.select()
            for button in self.link_inf_surface_bridges:
                button.select()
        if self.is_gln_selected.get():
            self.display_gln_content()

        if self.link_inf_is_smoothed.get():
            self.link_inf_smooth_button.select()

            bridges = self.selected_chains_links[filechain]
            hashcode = self.output_data[bridges][1]
            path_to_smooth_dir = self.link_directory + os.sep + hashcode + os.sep + "_sm" + os.sep
            ions_in_links = False
            closure = str(self.probabilistic_surface_idx[hashcode + str(self.displayed_link)]) \
                if self.notebook.getcurselection() == "Probabilistic" else "0"

            cmd.hide(representation="everything", selection="all")
            for i, elem in enumerate(bridges):
                atom = "ca"
                elem = elem.split(" ")
                if self.is_macrolink:
                    filename = elem[0]
                    sm_chain_xyz = hashcode + "_ch" + str(i) + "_" + filename + "_L"
                    color_filename = elem[0]
                    chain_in_file = "X"
                    chain_in_pymol = filename
                elif ".xyz" in elem[0]:
                    ions_in_links = True
                    filename = elem[0][:-4]
                    chain_in_file = elem[0].split(".pdb_")[1].split("_")[0]
                    chain_in_pymol = str(i)
                    sm_chain_xyz = hashcode + "_ch" + str(i) + "_" + filename + "_L"
                    color_filename = elem[0].split(".pdb_")[0]
                else:
                    filename = elem[0][:-1]
                    chain_in_file = elem[0][-1]
                    chain_in_pymol = chain_in_file
                    sm_chain_xyz = hashcode + "_ch" + str(i) + "_" + filename + ".pdb_" + chain_in_pymol + "_L"
                    color_filename = filename

                bridge = [elem[2], elem[1]] if int(elem[1]) > int(elem[2]) else [elem[1], elem[2]]
                file_atoms = [elem[1], elem[2]]

                if self.notebook.getcurselection() == "Probabilistic":
                    smooth_lvl = str(self.can_be_smoothed[filechain + str(self.displayed_link)][1])
                else:
                    smooth_lvl = str(self.can_be_smoothed[filechain][1])

                sm_chain_xyz += file_atoms[0] + "-" + file_atoms[1] + "_cl" + closure + "_sm" + smooth_lvl + ".xyz"
                sm_chain_pdb = sm_chain_xyz[:-4] + ".pdb"
                if ".xyz" in bridge:
                    sm_surface_file = hashcode + "_" + "s" + str(i) + "_" + sm_chain_xyz[:-4] + "_L" + file_atoms[0] + \
                                      "-" + file_atoms[1] + "_cl0.jm"
                else:
                    sm_surface_file = "s" + str(i) + "_" + sm_chain_xyz[:-4] + "_L" + file_atoms[0] + "-" + file_atoms[1] + \
                                      "_cl0.jm"
                if self.is_ions:
                    pymol_filename = "SM_" + str(i)
                else:
                    pymol_filename = "SM_" + chain_in_pymol
                self.convert_xyz_to_pdb(path_to_smooth_dir + sm_chain_xyz, path_to_smooth_dir, sm_chain_pdb,
                                        chain_in_file)
                self.get_cgo_coordinates(path_to_smooth_dir + sm_surface_file)
                if pymol_filename not in cmd.get_names() or ions_in_links:
                    cmd.load(filename=path_to_smooth_dir + sm_chain_pdb, object=pymol_filename)
                    if self.is_macrolink:
                        color = self.colors_of_sequence[i]
                        count_atoms = int(bridge[0]) + cmd.count_atoms(pymol_filename) + 1
                        beg = int(bridge[0])
                    elif ions_in_links:
                        color = util._color_cycle[self._chains[self._filenames.index(color_filename + ".pdb")].index(str(chain_in_file))]
                        count_atoms = self._marginal_atoms[filechain[:-1]][0] + cmd.count_atoms(pymol_filename) + 1
                    else:
                        color = util._color_cycle[self._chains[self._filenames.index(color_filename + ".pdb")].index(str(chain_in_pymol))]
                        count_atoms = self._marginal_atoms[filechain[:-1]][0] + cmd.count_atoms(pymol_filename) + 1
                        beg = self._marginal_atoms[filechain[:-1]][0] - 1
                    cmd.color(color=color, selection=pymol_filename)

                    if not self.is_macrolink:
                        if ions_in_links:
                            atoms = {'atoms': []}
                            cmd.iterate_state(state=1, selection="all and SM_" + str(i), expression="atoms.append(resi)",
                                              space=atoms)
                            straight_line = atoms['atoms']
                            for at1, at2 in zip(straight_line, straight_line[1:]):
                                cmd.bond(atom1=pymol_filename + " and id " + str(at1),
                                         atom2=pymol_filename + " and id " + str(at2))
                            cmd.bond(atom1=pymol_filename + " and id " + str(straight_line[-1]),
                                     atom2=pymol_filename + " and id " + str(straight_line[0]))
                        else:
                            straight_line = [beg, count_atoms] if beg < count_atoms else [count_atoms, beg]
                            for id in range(int(straight_line[0]), int(straight_line[1])):
                                cmd.bond(atom1=pymol_filename + " and id " + str(id),
                                         atom2=pymol_filename + " and id " + str(id + 1))

                if pymol_filename in cmd.get_names() and self.notebook.getcurselection() == "Probabilistic":
                    cmd.hide(representation="lines", selection="SM_" + chain_in_pymol)
                if "SEQ_SM_" + chain_in_pymol + "_*" in cmd.get_names():
                    cmd.show(representation="lines", selection="SEQ_SM_" + chain_in_pymol + "_*")

                if ions_in_links:
                    cmd.select("SEQ_SM_" + chain_in_pymol + "_" + bridge[0] + "_" + bridge[1],
                                selection="SM_" + chain_in_pymol)
                else:
                    cmd.select("SEQ_SM_" + chain_in_pymol + "_" + bridge[0] + "_" + bridge[1],
                                selection=pymol_filename + " and residue " + bridge[0] + "-" + bridge[1])
                cmd.show(representation="lines", selection="SEQ_SM_" + chain_in_pymol + "_" + bridge[0] + "_" + bridge[1])

                self.pymol_draw_piercings(i, "SM_" + chain_in_pymol + "_" + bridge[0] + "_" + bridge[1])
                self.pymol_draw_surface(i, "SM_" + chain_in_pymol + "_" + bridge[0] + "_" + bridge[1])

                if self.notebook.getcurselection() == "Probabilistic":
                    sm_br_cylinders = hashcode + "_ch" + str(i) + "ADD_" + filename + ".pdb_" + chain_in_file + "_L" \
                                      + bridge[0] + "-" + bridge[1] + "_cl" + closure + "_sm" + smooth_lvl + ".xyz"
                    points = open(path_to_smooth_dir + sm_br_cylinders).readlines()
                    self.br_cylinders = []
                    middle = points[1].rstrip().split(" ")[1:]
                    for point in [points[0], points[2]]:
                        point = point.rstrip().split(" ") + middle
                        line = []
                        for obj in point[1:]:
                            line.append(str("{0:.2f}".format(float(obj))))
                        self.br_cylinders.append(line)
                    self.pymol_draw_cylinders("SM_BR_" + chain_in_pymol + "_" + bridge[0] + "_" + bridge[1])
                else:
                    if not self.is_macrolink:
                        cmd.select("SM_BR_" + chain_in_pymol + "_" + bridge[0] + "_" + bridge[1],
                                   selection="(" + pymol_filename + " and id " + bridge[0] + " and name " + atom + ")+"
                                             + "(" + pymol_filename + " and id " + bridge[1] + " and name " + atom + ")")
                        cmd.bond(atom1=pymol_filename + " and id " + bridge[0],
                                 atom2=pymol_filename + " and id " + bridge[1])
                        cmd.show(representation='sphere', selection="SM_BR_" + chain_in_pymol + "_" + bridge[0] + "_" + bridge[1])
                        cmd.show(representation="sticks", selection="SM_BR_" + chain_in_pymol + "_" + bridge[0] + "_" + bridge[1])
                if not self.is_macrolink:
                    cmd.color(color=self.colors_of_sequence[i], selection="SEQ_SM_" + chain_in_pymol + "_" + bridge[0]
                                                                          + "_" + bridge[1])
                cmd.deselect()
            cmd.center(selection=pymol_filename)
            cmd.orient(selection=pymol_filename)
            cmd.clip(mode="slab", distance="1000")
            cmd.set(name="line_width", value="4")
            cmd.set(name="two_sided_lighting", value=1)
            cmd.set(name="sphere_color", value="orange", selection="all")
            cmd.set(name="stick_color", value="orange", selection="all")
            print "  Smoothed chain and bridges drawn in PyMOL..."
        else:
            cmd.delete(name="SM_*")
            cmd.delete(name="SURF_SM_*")
            cmd.delete(name="PIERC_SM_*")
            self.link_inf_smooth_button.deselect()
            if not self.is_macrolink:
                self.pymol_display_chain(self.selected_chains_links[self.displayed_filechain])
            self.pymol_display_link(self.displayed_link)

    def display_gln_content(self):
        if hasattr(self, "gln_matrices") and self.gln_matrices.winfo_exists():
            self.gln_matrices.grid_forget()
        self.gln_matrices = Pmw.Group(self.window_parent, tag_text='Gaussian Linking Number Matrices')
        self.gln_matrices.grid(column=0, row=2, pady=5, padx=5)

        rows = []
        for idx, i in enumerate(["Pair of components", "GLN Total", "GLN Average"]):
            rows.append(
                tk.Label(self.gln_matrices.interior(), text=i, width=self.gln_row_name_width, padx=0, pady=0,
                         font=self.row_name_font, justify="center"))
            rows[-1].grid(column=idx + 1, row=0)

        if self.is_macrolink or self._is_file_xyz[self.displayed_filechain[:-1]]:
            gln_filename = self._filenames[0][:-4]
        else:
            gln_filename = self.displayed_filechain[:-1]
        hashcode = self.output_data[self.selected_chains_links[self.displayed_filechain]][1]
        path_to_gln = self._full_path_to_dir + os.sep + gln_filename + os.sep + hashcode + os.sep + hashcode + \
                      "_glns.txt"

        self.gln_array_of_results = []
        height = 1
        content = open(path_to_gln, "r").readlines()
        glns, components = self.retrieve_gln_components_data(content)
        for elem in glns:
            row_element = []
            for obj in elem:
                row_element.append(
                    tk.Text(self.gln_matrices.interior(), bg="white", height=height + 1, padx=0, pady=0, bd=0,
                            width=self.gln_row_name_width+11, wrap="word", highlightthickness=0))
                try:
                    self.configure_tkinter_text(row_element[-1], float(obj), 1)
                except ValueError:
                    self.configure_tkinter_text(row_element[-1], obj, 1)
            row_element[0].configure(width=self.gln_row_name_width+10)
            self.gln_array_of_results.append(row_element)

        row = 1
        for filechain in self.gln_array_of_results:
            for idx, elem in enumerate(filechain):
                elem.grid(column=idx + 1, row=row)
            row += 1

        view_matrix_btns = []
        for idx in range(len(self.gln_array_of_results)):
            view_matrix_btns.append(tk.Button(self.gln_matrices.interior(), text="view matrix",
                                              justify="center", width=self.view_details_btn_width,
                                              command=lambda x=idx, y=components[idx]: self.display_gln_matrices(x, y)))
            view_matrix_btns[-1].grid(column=0, row=idx + 1)
        if len(view_matrix_btns) == 1:
            view_matrix_btns[0].invoke()

    def retrieve_gln_components_data(self, content):
        comp_dict = {}
        full_comp_data = {}
        components = content[0].rstrip().split(" ")[5:]
        gln_comp = []

        for idx, bridge in enumerate(components):
            bridge = bridge.split("_L")
            comp_dict[idx] = bridge[0][-1] + "" + bridge[1]
        for idx, line in enumerate(content[1:]):
            line = line.rstrip().split(" ")
            gln_comp.append([comp_dict[int(line[1])] + ", " + comp_dict[int(line[2])], line[4], line[-1]])
            full_comp_data[idx] = [components[int(line[1])], components[int(line[2])]]

        return gln_comp, full_comp_data

    def display_gln_matrices(self, x, content):
        self.displayed_gln = x
        self.color_line_in_array_of_results(self.gln_array_of_results, self.displayed_gln)

        hashcode = self.output_data[self.selected_chains_links[self.displayed_filechain]][1]
        if hasattr(self, "link_inf_is_smoothed") and self.link_inf_is_smoothed.get():
            target_dir = self.link_directory + os.sep + hashcode + os.sep + "_sm"
            working_dir = self.link_directory + os.sep + hashcode + os.sep + "_sm"
        else:
            target_dir = self.link_directory + os.sep + hashcode
            working_dir = self._full_path_to_dir if self.is_ions else self.link_directory
        if self.notebook.getcurselection() == "Probabilistic":
            command = self.generate_gln_probabilistic(content, working_dir)
        else:
            command = self.generate_gln_deterministic(content, working_dir)

        try:
            os.popen(command)
            if self.is_ions:
                self.separate_files_to_directory(self._full_path_to_dir, target_dir, hashcode + "_GLN_.*.png")
            self.separate_files_to_directory(self.current_working_dir, target_dir, hashcode + "_GLN_.*.png")
        except Exception:
            self.raise_popup_menu("The GLN value is equal to 0. No GLN files were generated.")

        if self.screen_height > 800:
            win_gln_matrices = Pmw.Group(self.gln_matrices.interior())
        else:
            win_gln_matrices = Pmw.ScrolledFrame(self.gln_matrices.interior(), usehullsize=1, hull_width=850, hull_height=300)
        win_gln_matrices.grid(column=0, row=len(self.gln_array_of_results)+2, columnspan=4)
        gln_matrices = mplt.figure.Figure(figsize=self.gln_fig_size, dpi=self.gln_fig_dpi, facecolor='white')

        files = command.split(" ")
        if self.is_macrolink:
            if hasattr(self, "link_inf_is_smoothed") and self.link_inf_is_smoothed.get():
                file1 = files[2].split("\\")[-1][:-4].split("macro_")[1].split("_L")[0]
                file2 = files[5].split("\\")[-1][:-4].split("macro_")[1].split("_L")[0]
            else:
                file1 = files[2].split("\\")[-1][:-4].split("macro_")[1]
                file2 = files[5].split("\\")[-1][:-4].split("macro_")[1]
        elif hasattr(self, "link_inf_is_smoothed") and self.link_inf_is_smoothed.get():
            file1 = files[2].split("\\")[-1].split(".pdb_")[1][0]
            if self.notebook.getcurselection() == "Probabilistic":
                file2 = files[3].split("\\")[-1].split(".pdb_")[1][0]
            else:
                file2 = files[5].split("\\")[-1].split(".pdb_")[1][0]
        else:
            if self.is_ions:
                file1 = files[2].split("\\")[-1][:-4].split(".pdb_")[1][0]
            else:
                file1 = files[2].split("\\")[-1][:-4][-1]
            file2 = files[5].split("\\")[-1][:-4].split(".pdb_")[1][0]
        name1 = hashcode + "_GLN_1" + file1 + "-2" + file2 + ".png"
        name2 = hashcode + "_GLN_2" + file2 + "-1" + file1 + ".png"

        for idx, img in enumerate([name1, name2]):
            self.fix_gln_matrices(gln_matrices, idx)
            img = mpimg.imread(target_dir + os.sep + img)
            gln_matrices.axes[idx].imshow(img, zorder=0)

        gln_matrices_canvas = FigureCanvasTkAgg(gln_matrices, master=win_gln_matrices.interior())
        gln_matrices_canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        gln_matrices_canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        gln_matrices_canvas.show()
        print "  GLN matrices displayed..."

    def generate_gln_deterministic(self, cnt, link_dir):
        comb = self.selected_chains_links[self.displayed_filechain]
        hashcode = self.output_data[self.selected_chains_links[self.displayed_filechain]][1]
        close = "11"
        args = ""
        ion = " -ion 1" if self.is_ions else ""

        for idx, components in enumerate(cnt):
            component = components.split("_L")
            br = self.fix_bridges(component[1])
            if self.is_ions:
                ch = str(idx)
            else:
                fixed_cmp = self.fix_component_name(components)
                ch = str(list(comb).index(fixed_cmp))
            if hasattr(self, "link_inf_is_smoothed") and self.link_inf_is_smoothed.get():
                args += link_dir + os.sep + hashcode + "_ch" + ch + "_" + components + "_cl0_sm" + \
                        str(self.smooth_level.get()) + ".xyz " + br[0] + " " + br[1] + " "
            else:
                args += link_dir + os.sep + component[0] + ".xyz " + br[0] + " " + br[1] + " "

        return self.prog_gln + args + "-out 0 -close " + close + ion + " -ht " + hashcode

    def fix_component_name(self, cmp):
        fixed_cmp = cmp.replace(".pdb_", "").replace("_L", " ")
        if fixed_cmp.count("-") > 1:
            fixed_cmp = fixed_cmp[::-1].replace("-", " ", 1)[::-1]
        else:
            fixed_cmp = fixed_cmp.replace("-", " ")
        return fixed_cmp

    def generate_gln_probabilistic(self, cnt, dir):
        comb = self.selected_chains_links[self.displayed_filechain]
        hashcode = self.output_data[comb][1]
        close = "00"
        args = ""

        for idx, components in enumerate(cnt):
            component = components.split("_L")
            br = self.fix_bridges(component[1])
            fixed_cmp = self.fix_component_name(components)
            ch = str(list(comb).index(fixed_cmp))
            if hasattr(self, "link_inf_is_smoothed") and self.link_inf_is_smoothed.get():
                cl = "0"
                args += dir + os.sep + hashcode + "_ch" + str(ch) + "_" + components + "_cl" + \
                        cl + "_sm" + str(self.smooth_level.get()) + ".xyz" + " "
            else:
                args += dir + os.sep + component[0] + ".xyz " + br[0] + " " + br[1] + " "

        return self.prog_gln + args + "-out 0 -close " + close + " -ht " + hashcode

    def fix_gln_matrices(self, gln, id):
        gln.add_axes([-0.28 + (id * 0.5), 0.4, 0.6, 0.3])
        gln.axes[id].set_position([-0.28 + (id * 0.5), 0, 1, 1])
        gln.axes[id].get_yaxis().set_visible(False)
        gln.axes[id].get_xaxis().set_visible(False)
        gln.axes[id].spines["top"].set_visible(False)
        gln.axes[id].spines["bottom"].set_visible(False)
        gln.axes[id].spines["left"].set_visible(False)
        gln.axes[id].spines["right"].set_visible(False)
        gln.axes[id].patch.set_visible(False)
        gln.axes[id].images = []


def _calculate_links(idx, command, out_data, sel_chains_links, complex_chains_links):
    link = command
    hashcode = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(4))
    link += " -ht " + hashcode
    homfly = call_homfly(link)
    if "ERROR" not in homfly:
        call_poly(hashcode)
        combination = global_combinations[idx]
        num_comp = str(len(combination))
        ncuc = call_ncuc(hashcode, num_comp, 1)
        process = filter(len, ncuc.split(" "))
        if process and process[0] != "\n":
            for id in range(len(process[3:])):
                process[id + 3] = process[id + 3].replace("X", " ")

            links = []
            for sele in process[2:-1]:
                elem = sele.split("(")
                l = elem[0] + "x" + elem[1][:-2]
                links.append(l)

            link_list = command.split(" ")
            num_comp = int(link_list[1])
            if int(link_list[2]) == 0: #Det
                if int(link_list[3]) == 1: #whole chains
                    zipped_list = link_list[4:4 + num_comp]
                else: #own loops
                    zipped_list = link_list[4:4 + (3 * (num_comp)):3]
            else: #Prob
                zipped_list = link_list[5:5 + (3 * (num_comp)):3]

            macroname_regex = re.compile("macro_(\d\w(_)*)*")
            for f in zipped_list:
                f = f.split(os.sep)[-1].split(".pdb_")
                filechain = f[0][:-4] if len(f) == 1 else f[0] + f[1][0]

                biggest_probability = 0
                for bridge in combination:
                    if macroname_regex.match(bridge.split(" ")[0]): # for macrolinks
                        bridge = bridge.split(" ")[0]
                    elif "_xyz2" in bridge:
                        filechain = f[0] + f[1][0]
                    elif ".xyz" in bridge: # for ions
                        file = bridge.split(".pdb_")[0]
                        chains = bridge.split(".pdb_")[1].split("_")[0]
                        bridge = file + chains
                    if filechain in bridge:
                        if links[0].split("x")[1]:
                            link_probability = links[0].split("x")[1][:-1]
                            if int(link_probability) > biggest_probability and biggest_probability != 100 \
                                    or filechain not in sel_chains_links:
                                biggest_probability = int(link_probability)
                                sel_chains_links[filechain] = combination
                                out_data[combination] = [links, hashcode]
                            if int(link_probability) > 0 and "COMPLEX" in links[0]:
                                complex_chains_links[combination] = links[1]
            return out_data, sel_chains_links, complex_chains_links


def call_homfly(command):
    try:
        return os.popen(command).read()
    except Exception:
        print "Something went wrong with executable file. Please make sure you changed access permission to " \
              "it (can be obtained by typing in console chmod a+x __homfly)."


def call_poly(hc):
    try:
        emcode = hc + "_EMCode.txt"
        os.popen(plugin_path + os.sep + "poly " + emcode)
    except Exception:
        print "Something went wrong with executable file. Please make sure you changed access permission to " \
              "it (can be obtained by typing in console chmod a+x __poly)."


def call_ncuc(hc, nc, inters=0):
    try:
        lmknot = hc + "_LMKNOT.txt" + (" -inters " + hc + "_inters.txt" if inters else "")
        return os.popen(plugin_path + os.sep + "ncucLinks " + lmknot + " -X 1 -comp " + nc).read()
    except Exception:
        print "Something went wrong with executable file. Please make sure you changed access permission to " \
              "it (can be obtained by typing in console chmod a+x __ncucLinks)."