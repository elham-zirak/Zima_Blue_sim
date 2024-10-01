########################################### Import Modules in Python ##########################################################

import numpy as np                          # for scientific computing
import pandas as pd                         # for data analysis and creating attenuation coefficient dictionaries        
import spekpy as sp                         # for getting poly-chromatic X-rays
import tifffile as tf                       # for writing the image as tiff file

############################################ Import Modules in Python for GUI #################################################

import tkinter as tk                        # for graphical user interface

""" acess classes from tkinter module, ttk: Tk themed widget set,
messagebox shows dialogues for errors
filedialog for showing the contrast value and save image
"""
from tkinter import ttk, messagebox, filedialog    
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg     # create canvas and embed plots in it
from matplotlib.figure import Figure        # for creating figure object

##################################### defining attenuation coefficient dictionaries ###########################################

# arguments: XmuDat data directory, seperator is space and there is no header
df = pd.read_csv('../Code/attenuation_coefficeint_data.txt', sep = '\s+', header = None)
# define headers for columns in dataframe by assigning names for columns
df.columns = ['Energy', 'tissue', 'water', 'bone', 'calcium', 'iodinated water']
# define a dictionary of density of input materials in g/cm^3 (based on XmuDat)
density = {'tissue' :1, 'water' : 1, 'bone' : 1.85, 'calcium' : 1.55, 'iodinated water': 1.01}

# converting columns of df (dataframe) to list for further use in dictionary, note that XmuDat gives mu/density,
# thereby it is needed to multiply it to density of each material:
mu_tissue = list(df['tissue']  * density['tissue'])
mu_water = list(df['water'] * density['water'])
mu_bone = list(df['bone'] * density ['bone'])
mu_calcium = list(df['calcium'] * density['calcium'])
mu_iodinated_water = list(df['iodinated water'] * density['iodinated water'])
energy = list(df['Energy'])

# defining dictionaries for linear attenuation coeffiecients in 1/cm
# in these dictionaries, key is i (energies), value is j (attenuation coefficients), zip function makes a pair of E and mu
mu_tissue_dict = {float(i): j for i,j in zip(energy, mu_tissue)}
mu_water_dict = {float(i): j for i,j in zip(energy, mu_water)}
mu_bone_dict = {float(i): j for i,j in zip(energy, mu_bone)}
mu_calcium_dict = {float(i): j for i,j in zip(energy, mu_calcium)}
mu_iodinated_water_dict = {float(i): j for i,j in zip(energy, mu_iodinated_water)}

# defining refrence dictionary for mu, key is compsition material and values are attenuation coefficients in 1/cm
all_mu = {'tissue': mu_tissue_dict, 'water': mu_water_dict, 'bone': mu_bone_dict, 'calcium': mu_calcium_dict, 'iodinated water' : mu_iodinated_water_dict}

################################################ pre-set values #############################################################

preset_values = {'kvp_ent': 100, 'targ_ent': 'W', 'theta_ent': 12, 'mas_ent': 100, 'filter_ent': 'Al', 'xfilter_ent': 3,
                 'obj_ent': 'tissue', 'det_ent': 'iodinated water', 'sdd_ent': 150, 'ssd_ent': 100, 'pixel_size': 0.005, 'FOV_ent': 2,
                 'scatter_ent': 0, 'det_radius': 0.5, 'xdet_ent': 0.5, 'xobj_ent': 5}

###################################### defining window and frame for graphical user interface ################################

root = tk.Tk()                                     # create root window
root.title('Poly Zima Blue')                       # root window title and dimension
root.geometry('850x900')                           # set geometry (width x height)for main window
frm = ttk.Frame(root, padding = 10)                # create a frame widget inside the root window and include a 10 pixel padding
frm.grid()                                         # grid the frame into columns and rows

######################################### creating labels for GUI ##################################################

# tk.Label(parent, text, options)
# For headers 'bg' argument sets the background color as an option.
# grid(coulmn, row) is used for geometry management.

# labels for X-ray tube
tk.Label(frm, text = 'X-ray tube inputs', bg = '#18BBF7').grid(column = 0, row = 0)
tk.Label(frm, text = 'Kilovoltage peak (kV)').grid(column = 0, row = 1)
tk.Label(frm, text = 'Anode target material').grid(column = 0, row = 2)
tk.Label(frm, text = 'The effective anode angle (°)').grid(column = 0, row = 3)
tk.Label(frm, text = 'Tube load (mAs)').grid(column = 0, row = 4)
tk.Label(frm, text = 'Filteration material').grid(column = 0, row = 5)
tk.Label(frm, text = 'Filteration thickness (mm)').grid(column = 0, row = 6)

# labels for Phantom and Geometrial Properties
tk.Label(frm, text = 'Phantom and Geometrial Properties', bg = '#3fc8f9').grid(column = 2, row = 0)
tk.Label(frm, text = 'Object composition').grid(column = 2, row = 1)
tk.Label(frm, text = 'Object thickness (cm)').grid(column = 2, row = 2)
tk.Label(frm, text = 'Detail composition').grid(column = 2, row = 3)
tk.Label(frm, text = 'Detail thickness (cm)').grid(column = 2, row = 4)
tk.Label(frm, text = 'Detail radius (cm)').grid(column = 2, row = 5)
tk.Label(frm, text = 'Source-Skin Distance or SSD (cm)').grid(column = 2, row = 6)
tk.Label(frm, text = 'Source-Detector Distance or SDD (cm)').grid(column = 2, row = 7)

# labels for image formation inputs
tk.Label(frm, text = 'Image Formation Properties', bg = '#18BBF7').grid(column = 4, row = 0)
tk.Label(frm, text = 'Detector pixel size (cm)').grid(column = 4, row = 1)
tk.Label(frm, text = 'Field of View (FOV) (cm)').grid(column = 4, row = 2)
tk.Label(frm, text = 'Scatter to Primary Ratio').grid(column = 4, row = 3)

################################# Creating entry fields ###################################

# tk.Entry(parent, width)
# insert(index, str), insert(0, preset_values['name of variable']) inserts the preset value for kVp at index 0
# grid(coulmn, row) is used for geometry management.

''' entry fields for X-ray tube '''

# kilovatge peak (kVp) entry
kvp_ent = tk.Entry(frm, width = 10)
kvp_ent.insert(0, str(preset_values['kvp_ent']))
kvp_ent.grid(column = 1, row = 1)

# anode angle entry
theta_ent = tk.Entry(frm, width=10)
theta_ent.insert(0, str(preset_values['theta_ent']))
theta_ent.grid(column = 1, row = 3)

# milliamperage (mAs) entry
mas_ent = tk.Entry(frm, width=10)
mas_ent.insert(0, str(preset_values['mas_ent']))
mas_ent.grid(column = 1, row = 4)

# filter material entry
xfilter_ent = tk.Entry(frm, width=10)
xfilter_ent.insert(0, str(preset_values['xfilter_ent']))
xfilter_ent.grid(column = 1, row = 6)

''' entry fields for Phantom and Geometrial Properties '''

# object thickness entry
xobj_ent = tk.Entry(frm, width=10)
xobj_ent.insert(0, str(preset_values['xobj_ent']))
xobj_ent.grid(column = 3, row = 2)

# detail thickness entry
xdet_ent = tk.Entry(frm, width=10)
xdet_ent.insert(0, str(preset_values['xdet_ent']))
xdet_ent.grid(column = 3, row = 4)

# detail radius entry
det_radius = tk.Entry(frm, width = 10)
det_radius.insert(0, str(preset_values['det_radius']))
det_radius.grid(column = 3, row = 5)

# Source-Skin Distance (SSD) entry (source-object)
ssd_ent = tk.Entry(frm, width = 10)
ssd_ent.insert(0, str(preset_values['ssd_ent']))
ssd_ent.grid(column = 3, row = 6)

# Source-Detector Distance (SDD) entry
sdd_ent = tk.Entry(frm, width = 10)
sdd_ent.insert(0, str(preset_values['sdd_ent']))
sdd_ent.grid(column = 3, row = 7)

''' entry fields for image formation '''

# pixel size entry
pixel_size = tk.Entry(frm, width=10)
pixel_size.insert(0, str(preset_values['pixel_size']))
pixel_size.grid(column = 5, row = 1)

# Field Of View (FOV) entry
FOV_ent = tk.Entry(frm, width=10)
FOV_ent.insert(0, str(preset_values['FOV_ent']))
FOV_ent.grid(column = 5, row = 2)

################################# Creating Option menus ####################################

# defining list of choices for option menus
target_lst = ['W', 'Mo', 'Rh']                                      # target material, restricted by SpekPy limitations 
filter_lst = ['Al', 'Mo', 'Sn', 'Cu']                               # filter material, 
object_lst = ['tissue']                                             # bject composition
detail_lst = ['water', 'bone', 'calcium', 'iodinated water']        # detail composition
R_lst = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]         # scatter to primary ratio list

# tk.StringVar(parent) returns a StringVar (control variable) instance for tk.OptionMenu() variable
# set() function 
# tk.OptionMenu (parent, variable, list of choices), the parent is the frame, variable shoudl be StrVariable instance and the list of choices are 

''' option menus for X-ray tube '''

# target material option menu
targ_ent = tk.StringVar(frm)
targ_ent.set(preset_values['targ_ent'])
targ_opt = tk.OptionMenu(frm, targ_ent, *target_lst)
targ_opt.grid(column = 1, row = 2)

# filter material option menu
filter_ent = tk.StringVar(frm)
filter_ent.set(preset_values['filter_ent'])
filter_opt = tk.OptionMenu(frm, filter_ent, *filter_lst)
filter_opt.grid(column = 1, row = 5)

''' option menus for Phantom and Geometrial Properties '''

# object material option menu
obj_ent = tk.StringVar(frm)
obj_ent.set(preset_values['obj_ent'])
obj_com = tk.OptionMenu(frm, obj_ent, *object_lst)
obj_com.grid(column = 3, row = 1)

# detail material option menu
det_ent = tk.StringVar(frm)
det_ent.set(preset_values['det_ent'])
det_com = tk.OptionMenu(frm, det_ent, *detail_lst)
det_com.grid(column = 3, row = 3)

''' option menu for image formation Properties '''

# scatter to primary ratio option menu
scatter_ent = tk.StringVar(frm)
scatter_ent.set(str(preset_values['scatter_ent']))
scatter_data = tk.OptionMenu(frm, scatter_ent, *R_lst)
scatter_data.grid(column = 5, row = 3)

############################################### claculate function  #######################################################

def calculate():
    
    # making the img_32bit, contrast, field and matrix size global variables to use them in other functions.
    global img_32bit, contrast, field_of_view, matrix

    # note that the entries from GUI are string variables,
    # so for the numerical variables it is needed to convert them to numerical datatypes such as float or int
    
    # entries for X-ray tube 
    kvp_data = float(kvp_ent.get())             # kilovoltage peak in kVp
    targ_data = targ_ent.get()                  # target material (can only be W, Mo, Rh based on spekpy limitations)
    theta = int(theta_ent.get())                # anode angle in degrees
    mas_data = int(mas_ent.get())               # milliamperage (mAs)
    filter_data = filter_ent.get()              # selected material for filter
    xfilter_data = float(xfilter_ent.get())     # thickness of filter in mm
    
    # entries for geometry and phantom properties
    sdd = float(sdd_ent.get())                  # source-detector distance
    ssd = float(ssd_ent.get())                  # source-skin (object) distance
    obj_material = obj_ent.get()                # object type
    det_material = det_ent.get()                # detail type
    det_rad_cm = float(det_radius.get())        # detail radius in cm
    xobj = float(xobj_ent.get())                # object thickness
    xdet = float(xdet_ent.get())                # detail thickness
    
    # previously we have defined all_mu dictionary.
    # now we define object and detail attenuation coefficient dictionaries based on the object and detail material selected by user in GUI.
    obj_dict = all_mu[obj_material]
    det_dict = all_mu[det_material] 
    
    # entries used for image formation
    pixel = float(pixel_size.get())                 # piexel size
    scatter_ratio = float(scatter_ent.get())        # scatter to primary ratio
    field_of_view = float(FOV_ent.get())            # field of view
    
    # basic calculations for image formation
    matrix = int(field_of_view/pixel)               # defining number of pixels in every direction of field of view
    det_rad_pixel = det_rad_cm/pixel                # converting detail radius from cm to pixel for image formation

    #################################################### Errors #############################################################
    
    # show error if kVp is not 10-100 kV for Mo target
    if float(kvp_data) not in range(10,101):
        messagebox.showerror('kVp Error', 'choose kilovoltage peak between 10 and 100 kV')
    
    # show error if kVp is not 20-50 kV for Mo target
    if float(kvp_data) not in range(20,51) and targ_data == 'Mo':
        messagebox.showerror('kVp Error', 'choose kilovoltage peak between 20 and 50 kV for Molybdenum target material')

    # show error if kVp is not 20-50 kV for Rh target
    if float(kvp_data) not in range(20,51) and targ_data == 'Rh':
        messagebox.showerror('kVp Error', 'choose kilovoltage peak between 20 and 50 kV for Rhodium target material')
    
    # show error if target angle is not in 1-180 range
    if theta not in range(7,21):
        messagebox.showerror('Anode angle Error', 'choose target angle between 7 and 20 degrees, with 1.0 increments')

    # setting condition that object thickness must be greater than detail thickness
    if xdet >= xobj:
        messagebox.showerror('Geometry Error', 'choose object thickness larger than detail thickness')
    
    # setting condition that SDD (source detector distance) must be greater than SSD (source skin distance)
    if ssd >= sdd:
        messagebox.showerror('Geometry Error', 'choose source skin distance smaller than source detector distance')
    
    ############################### getting X-ray spectra ####################################################
    
    # create X-ray spectrum with given kilovoltage peak, anode material and angle and tube load in distance z (SSD),
    # before entering the phantom with inherent filteration
    # Spek(kvp, th, targ, mas, x, y, z, physics, mu_data_source, dk)
    # A useful and handy resource on SpekPy: https://bitbucket.org/spekpy/spekpy_release/wiki/Home
    # filter('Al', 1) method accounts for inherent filtration of X-ray tube which is equivalent to 0.5-1 mm Al (see Farr's physics for medical imaging 1997- p29)
    s = sp.Spek(kvp = kvp_data, th = theta, targ = targ_data, mas = mas_data, z = ssd, physics = 'kqp', mu_data_source = 'pene', dk = 1).filter('Al', 1)
        
    # filteration: filter material and thickness are the arguments of filter method in SpekPy
    s.filter(matl = filter_data, t = xfilter_data)

    # get the spectrum, k is the energy and f is the fluence of photons (number of photons/ cm^2)
    k, f = s.get_spectrum(edges = True)

    # Using slicing method to remove duplicated values from k and f generated by get_spectrum method in SpekPy
    # [::2] syntax means [start:end:step], where start and end are omitted
    # We use this approach because the spectrum function may produce different numbers of values
    # depending on the X-ray tube input parameters.
    k = k[::2]
    f = f[::2]
                  
    # defining n_out with initial value of 0, which shows the number of photons outside of the detail
    n_out = 0                                            
    # defining n_in with initial value of 0, which shows the number of photons passing through the detail
    n_in = 0                                             

    ############################### fluence of X-ray spectrum (photons/cm^2) after attenuation ##################################
    
    # Calculate number of photons using Beer-Lambert's law
    # zip(k,f) means make a pair of calculated values for k and f in previous part,
    # where k is energy bin and f is correspodng fluences of X-ray spectrum
    # for loop iterates over k and f, i.e.
    # gives the energy bin (k) to obj_dict[i], which returns the value of linear attenuation coefficient based on energy and selected material by user
    # and gives the corresponding j to multply it by the exponentional part (np.exp : use numpy exponentional function for calculation)
    # n_out = number of photons per cm^2 outside the detail, which degrades exponentionally (Beer-Lambert's law)
    # n_in = number of photons per cm^2 passing through the detail
    
    for i, j in zip(k, f):
        n_out += j * np.exp(-obj_dict[i]*xobj)
        n_in += j * np.exp(-obj_dict[i]*(xobj-xdet)-(det_dict[i]*xdet))

 ######################################### Rescaling number of photons after attenuation ####################################

    # Rescaling using Inverse square law
    n_in = (ssd**2/sdd**2) * n_in
    n_out = (ssd**2/sdd**2) * n_out

    # rescaling to get number of photons per pixel
    n_in = (pixel**2) * n_in
    n_out = (pixel**2) * n_out
     
    # adding up scatter component    
    n_in = n_in + scatter_ratio *n_out
    n_out = n_out + scatter_ratio *n_out
    
    ########################################## creating coordinates; based on the matrix size #######################################
    
    # numpy arange retruns evenly spaced values within matrix size interval
    rrange = np.arange((-matrix/2),(matrix/2))               
    # using numpy tile function to divide the range into matrix number of tiles along first row to get x-coordinates
    x_coord = np.tile(rrange ,(matrix ,1))     
    # creating y coordinates by transposing x-coordinates
    y_coord = np.transpose(x_coord)
    
    # create img array using numpy where function; i.e:
    # place the total number of photons passing through the detail inside a circle with detail radius in units of pixel
    # and place the total number of photons in the background outside of the circle        
    img_out = np.where(x_coord**2 + y_coord**2 < det_rad_pixel**2, n_in, n_out)
    
    ################################################## add statistical noise ############################################# 
    
    # using np.random.poisson(lam) fo adding quantum mottle (noise) to the image due to the statistical fluctuations in the number of photons
    # lam (lambda) is indeed the expected number of events in a fixed-time interval. In X-ray imaging,
    # this corresponds to the expected number of photons hitting a particular detector pixel during the exposure time
    # lam = ii means setting the expected value in poisson function to the value of photons in each pixel by coordinates (x,y)
    img_ran = img_out                   # defining a new variable as img_ran which is initially the img_out
    for x in range(matrix):             # for every x in the range of matrix size
         for y in range(matrix):        # for every y in the range of matrix size (so we cover the whole image)        
            ii = img_out[x,y]           # get the original pixel value at coordinates (x,y) in the image 
            img_ran[x,y] = np.random.poisson(lam = ii)  # generate a random number from a Poisson distribution where lam is equal to ii and put it in img [x,y] (add to noise to whole pixels)

    ############################################### calculating contrast ##############################################
    
    # get the mean n_in where the coordinates are smaller than detail radius in pixel (inside the detail)
    n_in_avg = np.mean(img_ran[x_coord**2 + y_coord**2 < det_rad_pixel**2])
    
    # get the mean n_in where the coordinates are greater than or equal to detail radius in pixel (outside the detail)
    n_out_avg = np.mean(img_ran[x_coord**2 + y_coord**2 >= det_rad_pixel**2])
    
    # define local (Weber) contrast as the difference between average number of photons per pixel outside and inside the detail, divided by the background (outside the detail)
    contrast = ((n_out_avg - n_in_avg) / n_out_avg)
    
    # converting the imag_ran to numpy float32, i.e. a float of size 32 bits, which is desired for saving the image as tiff file   
    img_32bit = img_ran.astype(np.float32)
    
    return img_32bit, contrast, field_of_view, matrix       # return these variables when calling the main function via calculate button
  
################################################## Display image ########################################################

# create a function for activating when it is called by Display_image button to show the image inside a canvas

def Display_image():
    fig = Figure(figsize = (5, 5), dpi = 100)                   # create fig with a size of 5x5 inches and a resolution of 100 dots per inch.
    ax = fig.add_subplot(111)                                   # add a subplot to the figure. 111 means it's a 1x1 grid, and this is the first subplot.
    ax.imshow(img_32bit, cmap = 'gray')                         # display img_32bit on the subplot and set the colormap to grayscale 
    ax.axis('off')                                              # turn of axix labels
    canvas = FigureCanvasTkAgg(fig, master = frm)               # create a Tkinter canvas containing the matplotlib figure, and place it in frame (frm)   
    canvas_widget = canvas.get_tk_widget()                      # get the Tkinter widget corresponding to the canvas
    canvas_widget.grid(column = 0, row = 13, columnspan = 6)    # place the canvas widget in the Tkinter grid 0x13 with coulmnspan = 6
    canvas.draw()                                               # render img_32bit on the canvas 
    
############################################### showing contrast ###################################################

# define show_contrast function, so when the user clicks on contrast button, it calls this function.
def show_contrast():
    # set title and messge of the message box, use f string method to show txt and show the contrast value in 4 decimal places
    messagebox.showinfo(title= 'Contrast', message = f'The contrast is: {contrast:.4f}')

################################################ saving image ######################################################
def save_image():
    
    if 'img_32bit' not in globals():            # if the image is not generated previously
        messagebox.showerror('Image Error', 'No image to save. Please generate an image first.')    # show error message(title, txt)
        return                                  # exit this function if there is no image to save
    # open a file dialog for the user to choose where to save the image,
    # with default extension of the file as tiff, but the user can also save the file in any other file types
    file_path = filedialog.asksaveasfilename(defaultextension ='.tiff', filetypes = [('Tiff files', '*.tiff'), ('All files', '*.*')])
    
    if file_path:                               # if the filepath was selected by the user:
        tf.imwrite(file_path, img_32bit)        # save img_32bit in the filepath
        messagebox.showinfo("Success", f"Image saved successfully as 32-bit TIFF to {file_path}") # show this message after saving

############################################### define submit buttons ###############################################

# define Calculate button in the frame placed in 1*10 to activate calculate function
calculate_button = ttk.Button(frm, text = 'Calculate', command = calculate).grid(column = 1, row = 10)
# define image button in the frame placed in 2*10 to activate Display_image function
image_button = ttk.Button(frm, text = 'Display image', command = Display_image).grid(column = 2, row = 10)
# define contrast button in the frame placed in 3*10 to activate show_contrast function
contrast_button = ttk.Button(frm, text = 'Contrast', command = show_contrast).grid(column = 3, row = 10)
# define save image button in the frame placed in 4*10 to activate save_image function
save_button = ttk.Button(frm, text = 'Save Image', command = save_image).grid(column = 4, row = 10)
 
root.mainloop()         # keep the Tkinter window open and responsive to the actions done by user 
