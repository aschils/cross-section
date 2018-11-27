###############################
# Author: Arnaud Schils (NAPS)#
# arnaud.schils@gmail.com     #
###############################

import math
import numpy as np
import matplotlib.pyplot as plt
import sys
import Tkinter as tk
import tkFileDialog
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

AMU_TO_KG_FACTOR = 1.660539*10**-27
ELECTRON_CHARGE = 1.60217662*10**-19 #coulomb

main_window = tk.Tk()

class GuiEventHandler:

    def __init__(self):

        self.input_file_signal = ""
        self.input_file_noise = ""
        self.fig_id = 1

        self.nbr_of_events_with_noise = np.array([])
        self.noise = np.array([])
        self.I1 = 0
        self.I2 = 0
        self.q1 = 1
        self.q2 = -1
        self.m1 = 1
        self.m2 = 1
        self.E0_1 = 0
        self.E0_2 = 0
        self.V_obs_min = 0
        self.V_obs_max = 0
        self.eff_factor = 0
        self.form_factor = 0
        self.accept = 0
        self.delta_T = 0

        self.nbr_of_events = np.array([])
        self.V_obs = np.array([])
        self.delta_V_obs = 0
        self.E_rel = np.array([])
        self.v_rel_zeros_idx = np.array([])
        self.cross_section = np.array([])
        self.cross_section_error = np.array([])

        self.at_least_one_computation_done = False
        self.nbr_plot_windows = 0

    def get_nbr_plot_windows(self):
        return self.nbr_plot_windows
        
    #Properly close windows. If not called on window closing, error messages appear
    #when running the script on Windows 2000 OS
    def my_quit(self, window):
        window.destroy()
        self.nbr_plot_windows = self.nbr_plot_windows-1

    def valid_inputs(self):

        print("Signal file "+self.input_file_signal
        +" and noise file "+self.input_file_noise+"...")

        #Input file contains:
        #line_i represents an integer,
        #the number of events by the i_th canal during delta T seconds

        try:
            self.nbr_of_events_with_noise = np.loadtxt(self.input_file_signal)
        except:
            print("Invalid input file for signal.")
            return False
        try:
            self.noise = np.loadtxt(self.input_file_noise)
        except:
            print("Invalid input file for noise.")
            return False

        if self.nbr_of_events_with_noise.size != self.noise.size:
            print("Signal and noise arrays must have same sizes (check your input files).")
            return False

        #Currents I1 and I2 of beams 1 and 2 (nanoampere)
        try:
            self.I1 = float(entry_current1.get())*10**-9
        except ValueError:
            print("Invalid entry error: Beam 1 current must be numeric.")
            return False
        try:
            self.I2 = float(entry_current2.get())*10**-9
        except ValueError:
            print("Invalid entry error: Beam 2 current must be numeric.")
            return False
        #Ions charges
        try:
            self.q1 = float(entry_charge1.get())
        except ValueError:
            print("Invalid entry error: Beam 1 ion charge must be numeric.")
            return False
        try:
            self.q2 = float(entry_charge2.get())
        except ValueError:
            print("Invalid entry error: Beam 2 ion charge must be numeric.")
            return False
        charges_prod = self.q1*self.q2
        if charges_prod >= 0:
            print("Error: Ions must have opposite charges.")
            return False
        #Ions masses
        try:
            self.m1 = float(entry_mass1.get())
        except ValueError:
            print("Invalid entry error: Beam 1 ion mass must be numeric.")
            return False
        if(self.m1 <= 0):
            print("Invalid entry error: Beam 1 ion mass must be > 0.")
            return False
        try:
            self.m2 = float(entry_mass2.get())
        except ValueError:
            print("Invalid entry error: Beam 2 ion mass must be numeric.")
            return False
        if(self.m2 <= 0):
            print("Invalid entry error: Beam 2 ion mass must be > 0.")
            return False
        #Ions energies
        try:
            self.E0_1 = float(entry_energy1.get())
        except ValueError:
            print("Invalid entry error: Ion 1 energy must be numeric.")
            return False
        try:
            self.E0_2 = float(entry_energy2.get())
        except ValueError:
            print("Invalid entry error: Ion 2 energy must be numeric.")
            return False
        #Min and max observation potential in volt and absolute values
        try:
            self.V_obs_min = float(entry_vobs_min.get())
        except ValueError:
            print("Invalid entry error: Observation potential must be numeric.")
            return False
        try:
            self.V_obs_max = float(entry_vobs_max.get())
        except ValueError:
            print("Invalid entry error: Observation potential must be numeric.")
            return False
        if self.V_obs_max < self.V_obs_min:
            print("Max observation potential must be >= min observation potential.")
            return False
        try:
            self.eff_factor = float(entry_eff_factor.get())
        except ValueError:
            print("Invalid entry error: efficiency factor must be numeric.")
            return False
        if self.eff_factor <= 0:
            print("Efficiency must be > 0.")
            return False
        try:
            self.form_factor = float(entry_form_factor.get())
        except ValueError:
            print("Invalid entry error: form factor must be numeric.")
            return False
        if self.form_factor <= 0:
            print("Form factor must be > 0.")
            return False
        try:
            self.accept = float(entry_acceptance.get())
        except ValueError:
            print("Invalid entry error: acceptance must be numeric.")
            return False
        if self.accept <= 0:
            print("Acceptance must be > 0.")
            return False
        try:
            self.delta_T = float(entry_acquisition_time.get())
        except ValueError:
            print("Invalid entry error: acquisition time must be numeric.")
            return False
        if self.delta_T <= 0:
            print("Canal acquisition time must be > 0.")
            return False
        return True

    def compute_cross_section(self):

        self.nbr_of_events = self.nbr_of_events_with_noise-self.noise

        print("Number of events in canals:")
        print(self.nbr_of_events)

        #Error from Poisson law: sqrt(signal + 2*noise)
        abs_error_nbr_of_events_by_s = np.sqrt(self.nbr_of_events_with_noise
        +self.noise)/self.delta_T

        canal_idx = np.arange(0,self.nbr_of_events.size)
        nbr_of_canals = len(canal_idx)

        print("Number of canals: "+str(nbr_of_canals))

        if nbr_of_canals <= 0:
            print("At least one canal is required.")
            raise Exception("At least one canal is required.")

        #Convert canal number to observation potential V_obs
        if nbr_of_canals == 1:
            self.delta_V_obs = 0
        else:
            self.delta_V_obs = (self.V_obs_max-self.V_obs_min)/(nbr_of_canals-1.0)
        self.V_obs = self.V_obs_min+self.delta_V_obs*canal_idx

        print("Observation potentials Vobs (V):")
        print(self.V_obs)

        E1 = self.E0_1-self.q1*self.V_obs #eV
        E2 =  self.E0_2-self.q2*self.V_obs #eV

        print("Energies beam 1 (eV):")
        print(E1)
        print("Energies beam 2 (eV):")
        print(E2)

        max_nbr_of_events_idx = np.argmax(self.nbr_of_events)
        vobspeak = self.V_obs[max_nbr_of_events_idx]
        if max_nbr_of_events_idx > 0 and max_nbr_of_events_idx < (self.V_obs.size-1):
            #Correct energies: interpolation with parabola to find "true"
            #observation potential corresponding to max number of events
            x1 = self.V_obs[max_nbr_of_events_idx-1]
            x2 = self.V_obs[max_nbr_of_events_idx]
            x3 = self.V_obs[max_nbr_of_events_idx+1]
            y1 = self.nbr_of_events[max_nbr_of_events_idx-1]
            y2 = self.nbr_of_events[max_nbr_of_events_idx]
            y3 = self.nbr_of_events[max_nbr_of_events_idx+1]
            interpol_denominator = x1*(y3-y2)+x2*(y1-y3)+x3*(y2-y1)
            if interpol_denominator != 0:
                vobspeak = (x1**2*(y3-y2)+x2**2*(y1-y3)+x3**2*(y2-y1))/interpol_denominator/2.0

        E1_max_events = self.E0_1-self.q1*vobspeak
        E2_max_events = self.E0_2-self.q2*vobspeak
        de = (self.m1*E2_max_events-self.m2*E1_max_events)/(self.m1+self.m2)
        E1 = E1+de
        E2 = E2-de

        print("Energies beam 1 corrected (eV)")
        print(E1)
        print("Energies beam 2 corrected (eV)")
        print(E2)

        all_pos_energies = len(E1[E1 < 0]) == 0 and len(E2[E2 < 0]) == 0
        if not all_pos_energies:
            print("Error: negative energies forbidden.")
            raise Exception("Error: negative energies forbidden.")

        m1_kg = self.m1*AMU_TO_KG_FACTOR
        m2_kg = self.m2*AMU_TO_KG_FACTOR
        E1_joule = E1*ELECTRON_CHARGE
        E2_joule = E2*ELECTRON_CHARGE

        #E_cin = mv^2/2
        v1 = np.sqrt(2*E1_joule/m1_kg)*100.0 #cm/s
        v2 = np.sqrt(2*E2_joule/m2_kg)*100.0 #cm/s

        print("Speed beam 1 (cm/s)")
        print(v1)
        print("Speed beam 2 (cm/s)")
        print(v2)

        reduced_mass = self.m1*self.m2/(self.m1+self.m2)
        self.E_rel = reduced_mass*(E1/self.m1 + E2/self.m2 - 2*np.sqrt(E1*E2/(self.m1*self.m2)))

        print("Relative energies (eV):")
        print(self.E_rel)

        self.E_rel_joule = self.E_rel*ELECTRON_CHARGE
        reduced_mass_kg = reduced_mass*AMU_TO_KG_FACTOR
        v_rel = np.sqrt(2*np.abs(self.E_rel_joule)/reduced_mass_kg)*100.0 #cm/s

        print("Relative velocities (cm/s):")
        print(v_rel)

        nbr_of_events_by_s = self.nbr_of_events/self.delta_T

        #Remove 0 relative velocity points, avoid division by 0
        self.v_rel_zeros_idx = np.where(v_rel == 0)[0]
        v_rel = np.delete(v_rel, self.v_rel_zeros_idx)
        nbr_of_events_by_s = np.delete(nbr_of_events_by_s, self.v_rel_zeros_idx)
        v1 = np.delete(v1, self.v_rel_zeros_idx)
        v2 = np.delete(v2, self.v_rel_zeros_idx)
        abs_error_nbr_of_events_by_s = np.delete(abs_error_nbr_of_events_by_s,
        self.v_rel_zeros_idx)

        #!Points of infinite cross section not included in self.cross_section!
        #!Must be added afterwards in the draw plot methods!
        cross_section_div_nbr_of_events_by_s = self.q1*self.q2*ELECTRON_CHARGE**2*v1*v2/(self.I1*self.I2*v_rel*self.form_factor*self.eff_factor**2*self.accept)
        self.cross_section = np.abs(nbr_of_events_by_s*cross_section_div_nbr_of_events_by_s)

        print("Cross sections (cm^2):")
        print(self.cross_section)

        self.cross_section_error = np.abs(abs_error_nbr_of_events_by_s*cross_section_div_nbr_of_events_by_s)

        print("Error on cross sections (cm^2):")
        print(self.cross_section_error)

        self.at_least_one_computation_done = True

    def new_plot_fig(self):
        plot_fig = plt.figure(self.fig_id)
        self.fig_id = self.fig_id+1
        return plot_fig

    def create_plot_window(self, plot_fig):
        plt_window = tk.Tk()
        plt_window.protocol("WM_DELETE_WINDOW",
        lambda window=plt_window: self.my_quit(window))
        canvas = FigureCanvasTkAgg(plot_fig, master = plt_window)
        plot_widget = canvas.get_tk_widget()
        plot_widget.pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        toolbar = NavigationToolbar2Tk(canvas, plt_window)
        toolbar.update()
        canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        plot_fig.canvas.draw()
        self.nbr_plot_windows = self.nbr_plot_windows+1
        return plt_window

    def draw_plots(self):

        #Retrieve and check validity of parameters from graphical interface
        valid_inputs = self.valid_inputs()
        if not valid_inputs:
            print("Invalid user inputs. Cowardly giving up.")
            return
        try:
            self.compute_cross_section()
        except:
            print("Not able to compute cross section. Cowardly giving up.")
            return

        plot_fig_Vobs_cross_sec = self.new_plot_fig()
        V_obs_x = np.delete(self.V_obs, self.v_rel_zeros_idx)
        plt.errorbar(V_obs_x, self.cross_section, yerr=self.cross_section_error,
        ecolor="red")

        if V_obs_x.size > 0:
            plt.semilogy(V_obs_x, self.cross_section, "-4", color="blue")

        plt.xlabel("Observation potential [V]")
        plt.ylabel("Cross section [cm^2]")

        #Draw vertical red line at infinite cross section points
        for idx in self.v_rel_zeros_idx:
            #if idx <= max_nbr_of_events_idx:
            V_obs_i = self.V_obs_min+self.delta_V_obs*idx
            plt.axvline(V_obs_i, color="red")

        self.create_plot_window(plot_fig_Vobs_cross_sec)
        #self.window_plot1 = self.create_plot_window(plot_fig_Vobs_cross_sec)
        #self.window_plot1.protocol("wM_DELETE_WINDOW", self.quit_plot1)


        #After max value in nbr_of_events, data are non physical, remove them
        max_nbr_of_events_idx = np.argmax(self.nbr_of_events)
        idx_after_max = np.arange(max_nbr_of_events_idx,self.nbr_of_events.size)
        E_rel_x = np.delete(self.E_rel, idx_after_max)

        #For the cross section array, elements may already been removed due to
        #0 velocity. First remove from elements to remove already removed
        #elements.
        idx_to_remove_from_cross_sec = np.delete(idx_after_max,
        np.intersect1d(idx_after_max, self.v_rel_zeros_idx))
        #Then you need to shift remaining idx due to already removed elements
        idx_to_remove_from_cross_sec = idx_to_remove_from_cross_sec
        -self.v_rel_zeros_idx.size

        cross_section = np.delete(self.cross_section, idx_to_remove_from_cross_sec)
        cross_section_error = np.delete(self.cross_section_error,
            idx_to_remove_from_cross_sec)

        plot_fig_rel_energy_cross_sec = self.new_plot_fig()
        E_rel_x = np.delete(E_rel_x, self.v_rel_zeros_idx)
        plt.errorbar(E_rel_x , cross_section, yerr=cross_section_error, xerr=None,
        ecolor="red")
        if E_rel_x.size > 0:
            plt.loglog(E_rel_x , cross_section, "-4", color="blue")
        plt.xlabel("Relative energy [eV]")
        plt.ylabel("Cross section [cm^2]")

        #Draw vertical red line at infinite cross section points
        for idx in self.v_rel_zeros_idx:
            if idx <= max_nbr_of_events_idx:
                E_rel_i = self.E_rel[idx]
                plt.axvline(E_rel_i, color="red")

        self.create_plot_window(plot_fig_rel_energy_cross_sec)

    def browse_sig_file(self):
        selected_file_signal = tkFileDialog.askopenfilename(initialdir = "~",
        title = "select signal file")
        if selected_file_signal != "":
            self.input_file_signal = selected_file_signal

    def browse_noise_file(self):
        selected_file_noise = tkFileDialog.askopenfilename(initialdir = "~",
        title = "select noise file")
        if selected_file_noise != "":
            self.input_file_noise = selected_file_noise

    def write_to_file(self):

        if not self.at_least_one_computation_done:
            print("You must launch a computation before writing to file.")
            return

        try:
            out_file = tkFileDialog.asksaveasfile(initialdir = "~",
            title = "select output file", defaultextension=".txt")
        except:
            print("Unable to write to file.")
            return

        if out_file is None: #happens if dialog closed with "cancel"
            return

        # Vobs relative_energy nbr_of_events cross_section error_on_cross_section
        column_spacer = "    "
        out_file.write("Observation_potential_[V]"+column_spacer+
        "Relative_energy_[eV]"+column_spacer+"Number_of_events"+column_spacer+
        "Cross_section[cm^2]"+column_spacer+"Error_on_cross_section[cm^2]\n")

        nbr_of_0_v_rel_encountered = 0

        for i in range(0, self.V_obs.size):
            line = str(self.V_obs[i])+column_spacer+str(self.E_rel[i])
            line = line+column_spacer+str(self.nbr_of_events[i])+column_spacer

            if i in self.v_rel_zeros_idx:
                #then element for this observation potential has been removed
                #from cross section (infinite value...) and cross section error
                line = line+"infinity"+column_spacer+"NA"
                nbr_of_0_v_rel_encountered = nbr_of_0_v_rel_encountered+1
            else:
                corrected_idx = i-nbr_of_0_v_rel_encountered
                line = line+str(self.cross_section[corrected_idx])+column_spacer
                line = line+str(self.cross_section_error[corrected_idx])
            line = line+"\n"
            out_file.write(line)

        try:
            out_file.close()
        except:
            print("Unable to close file "+out_file)
            return
        print("Output file writed and closed.")

GEH = GuiEventHandler()

if len(sys.argv) == 3:
    GEH.input_file_signal = sys.argv[1]
    GEH.input_file_noise = sys.argv[2]

#Graphic interface Tkinter components
button_browse_sig_file = tk.Button(main_window, text="Browse signal file", command = GEH.browse_sig_file)
button_browse_sig_file.grid(row=0,column=1)

button_browse_noise_file = tk.Button(main_window, text="Browse noise file", command = GEH.browse_noise_file)
button_browse_noise_file.grid(row=0,column=4)

label_current1 = tk.Label(main_window, text="Beam 1 current [nA]:")
label_current1.grid(row=1,column=0)
entry_current1 = tk.Entry(main_window, bd =5)
entry_current1.grid(row=1,column=1)

label_current2 = tk.Label(main_window, text="Beam 2 current [nA]:")
label_current2.grid(row=2,column=0)
entry_current2 = tk.Entry(main_window, bd =5)
entry_current2.grid(row=2,column=1)

label_charge1 = tk.Label(main_window, text="Ion 1 charge [e]:")
label_charge1.grid(row=3,column=0)
entry_charge1 = tk.Entry(main_window, bd =5)
entry_charge1.grid(row=3,column=1)

label_charge2 = tk.Label(main_window, text="Ion 2 charge [e]:")
label_charge2.grid(row=4,column=0)
entry_charge2 = tk.Entry(main_window, bd =5)
entry_charge2.grid(row=4,column=1)

label_mass1 = tk.Label(main_window, text="Ion 1 mass [amu]:")
label_mass1.grid(row=5,column=0)
entry_mass1 = tk.Entry(main_window, bd =5)
entry_mass1.grid(row=5,column=1)

label_mass2 = tk.Label(main_window, text="Ion 2 mass [amu]:")
label_mass2.grid(row=6,column=0)
entry_mass2 = tk.Entry(main_window, bd =5)
entry_mass2.grid(row=6,column=1)

label_energy1 = tk.Label(main_window, text="Ion 1 energy [eV]:")
label_energy1.grid(row=7,column=0)
entry_energy1 = tk.Entry(main_window, bd =5)
entry_energy1.grid(row=7,column=1)

label_energy2 = tk.Label(main_window, text="Ion 2 energy [eV]:")
label_energy2.grid(row=1,column=3)
entry_energy2 = tk.Entry(main_window, bd =5)
entry_energy2.grid(row=1,column=4)

label_vobs_min = tk.Label(main_window, text="Observation potential min [V]:")
label_vobs_min.grid(row=2,column=3)
entry_vobs_min = tk.Entry(main_window, bd =5)
entry_vobs_min.grid(row=2,column=4)

label_vobs_max = tk.Label(main_window, text="Observation potential max [V]:")
label_vobs_max.grid(row=3,column=3)
entry_vobs_max = tk.Entry(main_window, bd =5)
entry_vobs_max.grid(row=3,column=4)

label_eff_factor = tk.Label(main_window, text="Detector efficiency")
label_eff_factor.grid(row=4,column=3)
entry_eff_factor = tk.Entry(main_window, bd =5)
entry_eff_factor.grid(row=4,column=4)

label_form_factor = tk.Label(main_window, text="Detector form factor")
label_form_factor.grid(row=5,column=3)
entry_form_factor = tk.Entry(main_window, bd =5)
entry_form_factor.grid(row=5,column=4)

label_acceptance = tk.Label(main_window, text="Detector acceptance")
label_acceptance.grid(row=6,column=3)
entry_acceptance = tk.Entry(main_window, bd =5)
entry_acceptance.grid(row=6,column=4)

label_acquisition_time = tk.Label(main_window, text="Detector acquisition time (per canal) [s]")
label_acquisition_time.grid(row=7,column=3)
entry_acquisition_time = tk.Entry(main_window, bd =5)
entry_acquisition_time.grid(row=7,column=4)

def go_button():
    global GEH
    GEH.draw_plots()
    global button_write_to_file
    button_write_to_file.config(state="normal")

button_go = tk.Button(main_window, text="Draw plots", command = go_button)
button_go.grid(row=14,column=2)

button_write_to_file = tk.Button(main_window, text="Write to file",
command=GEH.write_to_file, state="disabled")
button_write_to_file.grid(row=14,column=3)

def main_quit():
    global GEH
    if GEH.get_nbr_plot_windows() > 0:
        print("Close all plot windows first.")
        return
    main_window.destroy()
    sys.exit()

#main_window.protocol("WM_DELETE_WINDOW", lambda window=main_window: my_quit(window))
main_window.protocol("WM_DELETE_WINDOW", main_quit)
main_window.mainloop()
