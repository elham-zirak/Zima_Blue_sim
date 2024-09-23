########################################### Import Modules in Python ###########################################################################

import numpy as np                      #import numpy for scientific computing
#importing math module for using mathematical functions defined by C standard
import math                             #import math for mathematical functions
#importing pyplot interface as plt for adjusting figures like matlab in matplotlib
from matplotlib import pyplot as plt

#from matplotlib import image as img             

#importing pandas module as pd for data analysis
import pandas as pd
#importing spekpy as sp for getting poly-chromatic x-rays                             
import spekpy as sp

import scipy as sc

#importing tkinter module for graphical user interface
import tkinter
#importing everything from tkinter module, so there is no need to prefix tkinter
from tkinter import *
#importing ttk (themed tk) from tkinter module, to access Tk themed widget set 
from tkinter import ttk
#importing messagebox from tkinter to show dialogues for errors
from tkinter import messagebox
#importing FigureCanvasTkAgg class to create canvas and embed plots in it, backends are like plugs between tkinter and matplotlib                                   
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
#importing Figure class from matplotlib to create figure object
from matplotlib.figure import Figure 
#importing fieldalog for showing the contrast value 
from tkinter import filedialog

##################################### defining attenuation coefficient dictionaries  ###############################################
#create a dataframe from csv file using pandas module
#arguments: XmuDat data directory, seperator is space and there is no header
df = pd.read_csv('../mono/attenuation_coefficeint_data.txt', sep='\s+', header = None) #arguments: XmuDat data directory, seperator is space and there is no header        
#define headers for columns in dataframe by assigning names for columns
df.columns = ['Energy', 'tissue', 'water', 'bone', 'calcium', 'Iodine water']

#define a dictionary of density of input materials in g/cm^3 (based on XmuDat)
density = {'tissue' :1, 'water' : 1, 'bone' : 1.85, 'calcium' : 1.55, 'Iodine water': 1.01}

#converting columns of df (dataframe) to list, note that XmuDat gives mu/density,
# thereby it is needed to multiply it to density of each material:
mu_tissue = list(df['tissue']  * density['tissue'])
mu_water = list(df['water'] * density['water'])
mu_bone = list(df['bone'] * density ['bone'])
mu_calcium = list(df['calcium'] * density['calcium'])
mu_iodine = list(df['Iodine water'] * density['Iodine water'])
energy = list(df['Energy'])

#defining dictionaries for linear attenuation coeffiecients in 1/cm
#in these dictionaries, key is i (energies), value is j (attenuation coefficients), zip function makes a pair of E and mu
mu_tissue_dict = {i: j for i,j in zip(energy, mu_tissue)}
mu_water_dict = {i: j for i,j in zip(energy, mu_water)}
mu_bone_dict = {i: j for i,j in zip(energy, mu_bone)}
mu_calcium_dict = {i: j for i,j in zip(energy, mu_calcium)}
mu_iodine_dict = {i: j for i,j in zip(energy, mu_iodine)}

#defining refrence dictionary for mu, key is compsition material and values are attenuation coefficients in 1/cm
all_mu = {'tissue': mu_tissue_dict, 'water': mu_water_dict, 'bone': mu_bone_dict, 'calcium': mu_calcium_dict, 'Iodine water' : mu_iodine_dict}

################################################ pre-set values ####################################################################

preset_values = {'kvp_ent': 40, 'targ_ent': 'W', 'theta_ent': 12, 'mas_ent': 100, 'filter_ent': 'Al', 'xfilter_ent': 3,
                 'obj_ent': 'tissue', 'det_ent': 'calcium', 'sdd_ent': 150, 'ssd_ent': 100, 'pixel_size': 0.005, 'FOV_ent': 2,
                 'scatter_ent': 0.1, 'det_radius': 0.5, 'xdet_ent': 0.5, 'xobj_ent': 5}

############################################# defining window for graphical user interface ###############################

root = Tk()                                                                     # create root window
root.title('Zima Blue')                                                         # root window title and dimension
root.geometry('800x700')                                                        # Set geometry (width x height)
frm = ttk.Frame(root, padding=10)                                               # create a frame inside the root window with padding 10
frm.grid()                                                                      # grid the frame into columns and rows

############################## defining labels for x-ray tube ##################################################

Label(frm, text = 'x-ray tube inputs', bg = 'cornsilk').grid(column = 0, row = 0)       #header in 0*0
Label(frm, text = 'kilovoltage peak (kV)').grid(column = 0, row = 1)                    #kvp label in 0*1
Label(frm, text = 'anode target material').grid(column = 0, row = 2)                    #target material label in 0*2
Label(frm, text = 'The effective anode angle (°)').grid(column = 0, row = 3)            #target angle label in 0*3
Label(frm, text = 'Tube load (mAs)').grid(column = 0, row = 4)                          #tube load label in 0*4
Label(frm, text = 'Filteration material').grid(column = 0, row = 5)                     #filter material label in 0*5
Label(frm, text = 'Filteration thickness (mm)').grid(column = 0, row = 6)               #filter thickness label in 0*6

############################## defining option menus and entry boxes for x-ray tube ###########################

target_lst = ['W', 'Mo', 'Rh']                                  #create a list for use in option menu of otarget material
filter_lst = ['Al', 'Mo', 'Sn', 'Cu']                           #create a list for use in option menu of filter material

# create entry for kvp in 1*1
kvp_ent = Entry(frm, width=10)
kvp_ent.insert(0, str(preset_values['kvp_ent']))
kvp_ent.grid(column = 1, row = 1)

# create option menu for target material in 1*2
targ_ent = tkinter.StringVar(frm)                                #create a string variable for target material
targ_ent.set(preset_values['targ_ent'])
targ_opt = OptionMenu(frm, targ_ent, *target_lst)               #create an option menu & assign str variables in target_lst to it
targ_opt.grid(column = 1, row = 2)                               #put det com options in 3*2

# create entry for target angle in 1*3
theta_ent = Entry(frm, width=10)
theta_ent.insert(0, str(preset_values['theta_ent']))
theta_ent.grid(column = 1, row = 3)

# create entry for tube load in 1*4
mas_ent = Entry(frm, width=10)
mas_ent.insert(0, str(preset_values['mas_ent']))
mas_ent.grid(column = 1, row = 4)

# create option menu for filter material in 1*5
filter_ent = tkinter.StringVar(frm)                           #create a string variable for targ material
filter_ent.set(preset_values['filter_ent'])
filter_opt = OptionMenu(frm, filter_ent, *filter_lst)         #create an option menu & assign str variables in filter_lst to it
filter_opt.grid(column = 1, row = 5)                               #put filter options in 1*5

# create entry for tube load in 1*6
xfilter_ent = Entry(frm, width=10)
xfilter_ent.insert(0, str(preset_values['xfilter_ent']))
xfilter_ent.grid(column = 1, row = 6)

############################## defining labels for geometries ########################################################

Label(frm, text = 'Geometry', bg = 'mistyrose').grid(column = 2, row = 0)         #header in 2*0
Label(frm, text = 'object composition').grid(column = 2, row = 1)                   #obj com label in 2*1
Label(frm, text = 'object thickness (cm)').grid(column = 2, row = 2)               #xobj label in 2*2
Label(frm, text = 'detail composition').grid(column = 2, row = 3)                   #det com label in 2*3
Label(frm, text = 'detail thickness (cm)').grid(column = 2, row = 4)               #xdet label in 2*4
Label(frm, text = 'detail radius (cm)').grid(column = 2, row = 5)                  #det radius label in 2*5
Label(frm, text = 'Source to Skin Distance or SSD (cm)').grid(column = 2, row = 6)                 #SSD label in 2*6
Label(frm, text = 'Source to Detector Distance or SDD (cm)').grid(column = 2, row = 7)       #SDD label in 2*7

###################### defining option menus and entry boxes for object and detail in column 3 ###################################
object_lst = ['tissue']           #create a list for use in option menus of object
detail_lst = ['water', 'bone', 'calcium', 'Iodine water']           #create a list for use in option menus of object

# create option menu for object composition
obj_ent = tkinter.StringVar(frm)                           #create a string variable for obj com
obj_ent.set(preset_values['obj_ent'])
obj_com = OptionMenu(frm, obj_ent, *object_lst)          #create an option menu & assig str variables in material_lst to it 
obj_com.grid(column = 3, row = 1)                               #put obj com options in 3*1

# create entry for object thickness in 3*2
xobj_ent = Entry(frm, width=10)
xobj_ent.insert(0, str(preset_values['xobj_ent']))
xobj_ent.grid(column = 3, row = 2)

# create option menu for detail composition
det_ent = tkinter.StringVar(frm)                           #create a string variable for det com
det_ent.set(preset_values['det_ent'])
det_com = OptionMenu(frm, det_ent, *detail_lst)          #create an option menu & assign str variables in material_lst to it
det_com.grid(column = 3, row = 3)                               #put det com options in 3*3

# create entry for object thickness in 3*4
xdet_ent = Entry(frm, width=10)
xdet_ent.insert(0, str(preset_values['xdet_ent']))
xdet_ent.grid(column = 3, row = 4)

# create entry for detail radius in 3*5
det_radius = Entry(frm, width=10)
det_radius.insert(0, str(preset_values['det_radius']))
det_radius.grid(column = 3, row = 5)

# create entry for SSD (source skin distance) in 3*6
ssd_ent = Entry(frm, width=10)
ssd_ent.insert(0, str(preset_values['ssd_ent']))
ssd_ent.grid(column = 3, row = 6)

# create entry for SDD (source detector distance) in 3*7
sdd_ent = Entry(frm, width=10)
sdd_ent.insert(0, str(preset_values['sdd_ent']))
sdd_ent.grid(column = 3, row = 7)

############################## defining labels for image ########################################################

Label(frm, text = 'image properties', bg = 'lightblue').grid(column = 4, row = 0)                #header in 4*0
Label(frm, text = 'detector pixel size (cm)').grid(column = 4, row = 1)                     #pixel size label in 4*1
Label(frm, text = 'Field of View (FOV) (cm)').grid(column = 4, row = 2)                     #FOV label in 4*2
Label(frm, text = 'Scatter to Primary Ratio').grid(column = 4, row = 3)                     #scatter to primary ratio label in 4*3

############################## defining option menus and entry boxes for image ###########################

scatter_lst = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]           #create a list for use in option menus of obj and det composition

# create entry for pixel size in 5*1
pixel_size = Entry(frm, width=10)
pixel_size.insert(0, str(preset_values['pixel_size']))
pixel_size.grid(column = 5, row = 1)

# create entry for FOV in 5*2
FOV_ent = Entry(frm, width=10)
FOV_ent.insert(0, str(preset_values['FOV_ent']))
FOV_ent.grid(column = 5, row = 2)

# create option menu for scatter to primary ratio
scatter_ent = tkinter.StringVar(frm)                                #create a string variable for scatter_ent
scatter_ent.set(str(preset_values['scatter_ent']))
scatter_data = OptionMenu(frm, scatter_ent, *scatter_lst)          #create an option menu & assign str variables in scatter_lst to it
scatter_data.grid(column = 5, row = 3)                              #put scatter to primary options in 5*3

############################################### main function  #######################################################

def calculate():
    # defining global variables
    global img_ran, contrast, canvas, ax, rs

    #entries for x-ray tube
    kvp_data = float(kvp_ent.get())
    targ_data = targ_ent.get()
    theta = int(theta_ent.get())
    mas_data = int(mas_ent.get())
    filter_data = filter_ent.get()
    xfilter_data = float(xfilter_ent.get())
    
    #entries for geometry
    sdd = float(sdd_ent.get())              #source to detector distance
    ssd = float(ssd_ent.get())              #source to skin (object) distance
    obj_material = obj_ent.get()            #object type
    det_material = det_ent.get()            #detail type 
    det_rad_cm = float(det_radius.get())        #detail radius
    xobj = float(xobj_ent.get())                #object thickness
    xdet = float(xdet_ent.get())                #detail thickness
    
    #define mu dictionaries for object and detail
    obj_dict = all_mu[obj_material]
    det_dict = all_mu[det_material] 
    
    #entries for image
    pixel = float(pixel_size.get())                 #getting pixel size from entries
    scatter_ratio = float(scatter_ent.get())              #getting scatter to primary ratio for calculations
    field = float(FOV_ent.get())                    #defining field of view
    
    #defining global variables for contrast and final image
    global contrast
    global img_ran
    global canvas

    #calculation
    matrix = int(field/pixel)                       #defining number of pixels in every direction of field of view
    det_rad_pixel = det_rad_cm/pixel                #converting detail radius from cm to pixel for image formation

    
    ############################### getting x-ray spectra ####################################################
    #create xray spectrum with given kilovoltage peak, target material and angle and tube load in distance z (SSD),
    #before entering the phantom with inherent filteration
    s = sp.Spek(kvp = kvp_data, th = theta, targ = targ_data, mas = mas_data, z = ssd, mu_data_source = 'nist', dk = 1)
        
    #filteration: filter material and thickness
    s.filter(matl = filter_data, t = xfilter_data)

    #get the spectrum, k is the energy and f is the fluence of photons in SSD
    k, f = s.get_spectrum(edges=True)

    k = k[::2]                                          #removing duplicate values from k and converting ndarray to list
    f = f[::2]                                          #removing duplicate values from f

    ###################################### Rescale constant ############################################
    a = (ssd**2/sdd**2) * (pixel**2)                    #define rescale constant so we get number of photons                    

    #setting initial value of n_in and n_out
    n_out = 0                                           
    n_in = 0                                            

    ####################################### get number of photons after attenuation ##################################
    #getting the number of photons when they arrive to detector 
    for i, j in zip(k, f):
        n_out += (a *j) * math.exp(-obj_dict[i]*xobj)
        n_in += (a *j) * math.exp(-obj_dict[i]*(xobj-xdet)-(det_dict[i]*xdet))

    ########################################### Scattered photons calculation #############################################
    # adding up scatter component    
    n_in_total = n_in + scatter_ratio *n_out                        #calculating total number of photons passing through the detail
    n_out_total = n_out + scatter_ratio *n_out                      #calculationg total number of photons in the background
    
    #################################################### Errors #############################################################

    #show error if kVp is not 20-50 kV for Mo target
    if float(kvp_data) not in range(20,51) and targ_data == 'Mo':
        messagebox.showerror('kVp Error', 'choose kilovoltage peak between 20 and 50 kV for Molybdenum target material')

    #show error if kVp is not 20-50 kV for Rh target
    if float(kvp_data) not in range(20,51) and targ_data == 'Rh':
        messagebox.showerror('kVp Error', 'choose kilovoltage peak between 20 and 50 kV for Rhodium target material')
    
    #show error if target angle is not in 1-180 range
    if theta not in range(7,21):
        messagebox.showerror('Anode angle Error', 'choose target angle between 7 and 20 degrees, with 1.0 increments')

    #setting condition that object thickness must be greater than detail thickness
    if xdet >= xobj:
        messagebox.showerror('Geometry Error', 'choose object thickness larger than detail thickness')
    
    #setting condition that SDD (source detector distance) must be greater than SSD (source skin distance)
    if ssd >= sdd:
        messagebox.showerror('Geometry Error', 'choose source skin distance smaller than source detector distance')
    
      
    ########################################### creating coordinates; based on the matrix size #######################################
    
    #numpy arange retruns evenly spaced values within matrix size interval
    rrange = np.arange((-matrix/2),(matrix/2))               
    #using numpy tile function to divide the range into matrix number of tiles along first row to get x-coordinates
    x_coord = np.tile(rrange ,(matrix ,1))     
    #creating y coordinates by transposing x-coordinates
    y_coord = np.transpose(x_coord)
    
    #create img array using numpy where function; i.e:
    #place the total number of photons passing through the detail inside a circle with detail radius in units of pixel
    # and place the total number of photons in the background outside of the circle        
    img_out = np.where(x_coord**2 + y_coord**2 < det_rad_pixel**2, int(n_in_total),int(n_out_total))
       
    # add statistical noise 
    img_ran = img_out
    for x in range(matrix):
         for y in range(matrix):
            ii = img_out[x,y]               #getting the original pixel value
            img_ran[x,y] = np.random.poisson(lam=ii)            #assigning 

    # Create a new Figure and add a subplot
    fig = Figure(figsize=(5, 4), dpi=100)
    ax = fig.add_subplot(111)

    # Display the image on the subplot
    ax.imshow(img_ran)
    ax.axis('off')

    # Create a FigureCanvasTkAgg object
    canvas = FigureCanvasTkAgg(fig, master = root)
    canvas_widget = canvas.get_tk_widget()
    canvas_widget.grid(column=0, row=11, columnspan=6)

    # Draw the canvas
    canvas.draw()

############################################### calculating contrast ##############################################
    n_in_avg = np.mean(img_ran[x_coord**2 + y_coord**2 < det_rad_pixel**2])
    n_out_avg = np.mean(img_ran[x_coord**2 + y_coord**2 >= det_rad_pixel**2])
    contrast = ((n_out_avg - n_in_avg) / n_out_avg)

############################################### showing contrast ###################################################
# def show_contrast():
    contrast_value = contrast  #contrast value
    messagebox.showinfo(title= 'Contrast', message = f'The contrast is: {contrast_value:.4f}')

################################################ saving image ######################################################
def save_image():
    if 'img_ran' not in globals():
        messagebox.showerror("Error", "No image to save. Please generate an image first.")
        return
    
    file_path = filedialog.asksaveasfilename(defaultextension=".tiff",
                                             filetypes=[("Tiff files", "*.tiff"),
                                                        ("All files", "*.*")])
    if file_path:
        plt.imsave(file_path, img_ran)
        messagebox.showinfo("Success", f"Image saved successfully to {file_path}")


############################################### define submit buttons ###############################################
submit_button_1 = ttk.Button(frm, text='Display image', command = calculate).grid(column = 1, row = 10)
contrast_button = ttk.Button(frm, text='Contrast', command = show_contrast).grid(column=3, row=10)
save_button = ttk.Button(frm, text='Save Image', command = save_image).grid(column=5, row=10)
 
root.mainloop()