########################################### Import Modules in Python ##########################################################

import numpy as np                          # for scientific computing
import pandas as pd                         # for data analysis and creating attenuation coefficient dictionaries        
import tifffile as tf                       # for writing the image as tiff file

############################################ Import Modules in Python for GUI #################################################

import tkinter as tk                        # for graphical user interface

""" acess classes from tkinter module, ttk: Tk themed widget set,
messagebox shows dialogues for errors
filedialog for showing the contrast value and save image
"""
from tkinter import ttk, messagebox, filedialog
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg     # create canvas and embed plots in it, used for saving image
from matplotlib.figure import Figure                                # for creating figure object used for saving image

##################################### defining attenuation coefficient dictionaries ###########################################

# data retrived from XmuDat, for tissue, water (representing a soft tomour), bone and calcium (representing calcification)
# and iodinated water (99% water + 1% iodine) representing a contrast agent
# arguments: XmuDat data directory, seperator is space and there is no header
df = pd.read_csv('attenuation_coefficeint_data.txt', sep='\s+', header = None)
# define headers for columns in dataframe by assigning names for columns
df.columns = ['Energy', 'tissue', 'water', 'bone', 'calcium', 'iodinated water']
# define a dictionary of density of input materials in g/cm³ (based on XmuDat)
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
mu_tissue_dict = {i: j for i,j in zip(energy, mu_tissue)}
mu_water_dict = {i: j for i,j in zip(energy, mu_water)}
mu_bone_dict = {i: j for i,j in zip(energy, mu_bone)}
mu_calcium_dict = {i: j for i,j in zip(energy, mu_calcium)}
mu_iodinated_dict = {i: j for i,j in zip(energy, mu_iodinated_water)}

# defining refrence dictionary for mu, key is compsition material and values are attenuation coefficients in 1/cm
all_mu = {'tissue': mu_tissue_dict, 'water': mu_water_dict, 'bone': mu_bone_dict, 'calcium': mu_calcium_dict, 'iodinated water' : mu_iodinated_dict}

################################################ pre-set values #############################################################

# defining a dictionary of pre-set values, dictionary = {'key_1': value_1, 'key_2': value_2, ....}
preset_values = { 'n0_ent': 100000000, 'mono_energy_ent': 20 ,'obj_ent': 'tissue', 'det_ent': 'calcium', 'sdd_ent': 60,
                 'ssd_ent': 50, 'pixel_size': 0.01, 'FOV_ent': 2, 'scatter_ent': 0, 'det_radius': 0.25, 'xdet_ent': 0.1,
                  'xobj_ent': 5}

###################################### defining window and frame for graphical user interface ################################

root = tk.Tk()                                                    # create root window
root.title('Mono Zima Blue')                                      # root window title and dimension
root.geometry('800x700')                                          # Set geometry (width x height)
frm = ttk.Frame(root, padding = 10)                               # create a frame inside the root window with padding 10
frm.grid()                                                        # grid the frame into columns and rows

######################################### creating labels for GUI ##################################################

# tk.Label(parent, text, options)
# For headers 'bg' argument sets the background color as an option.
# grid(coulmn, row) is used for geometry management.

# labels for X-ray 
tk.Label(frm, text = 'X-ray inputs', bg = '#18BBF7').grid(column = 0, row = 0)
tk.Label(frm, text = 'Energy (keV)').grid(column = 0, row = 1)
tk.Label(frm, text = 'Number of incident \nphotons/cm\u00B2 at source').grid(column = 0, row = 2)

# labels for Phantom and Geometrical Properties
tk.Label(frm, text = 'Phantom and Geometrical Propertie', bg = '#18BBF7').grid(column = 2, row = 0)
tk.Label(frm, text = 'object composition').grid(column = 2, row = 1)
tk.Label(frm, text = 'object thickness (cm)').grid(column = 2, row = 2)
tk.Label(frm, text = 'detail composition').grid(column = 2, row = 3)
tk.Label(frm, text = 'detail thickness (cm)').grid(column = 2, row = 4)
tk.Label(frm, text = 'detail radius (cm)').grid(column = 2, row = 5)
tk.Label(frm, text = 'Source to Skin Distance or SSD (cm)').grid(column = 2, row = 6)
tk.Label(frm, text = 'Source to Detector Distance or SDD (cm)').grid(column = 2, row = 7)

# labels for image formation inputs
tk.Label(frm, text = 'Image Formation Properties', bg = '#18BBF7').grid(column = 4, row = 0)
tk.Label(frm, text = 'detector pixel size (cm)').grid(column = 4, row = 1)
tk.Label(frm, text = 'Field of View (FOV) (cm)').grid(column = 4, row = 2)
tk.Label(frm, text = 'Scatter to Primary Ratio').grid(column = 4, row = 3)

################################# Creating entry fields ###################################

# tk.Entry(parent, width)
# insert(index, str), insert(0, preset_values['name of variable]) inserts the preset value for kVp at index 0
# grid(coulmn, row) is used for geometry management.

''' entry fields for X-ray tube '''

# energy (keV) entry
mono_energy_ent = tk.Entry(frm, width=10)
mono_energy_ent.insert(0, str(preset_values['mono_energy_ent']))
mono_energy_ent.grid(column = 1, row = 1)

# number of incident photons per cm² entry
n0_ent = tk.Entry(frm, width=10)
n0_ent.insert(0, str(preset_values['n0_ent']))
n0_ent.grid(column = 1, row = 2)

''' entry fields for Phantom and Geometrical Properties '''

# object thickness entry
xobj_ent = tk.Entry(frm, width=10)
xobj_ent.insert(0, str(preset_values['xobj_ent']))
xobj_ent.grid(column = 3, row = 2)

# detail thickness entry
xdet_ent = tk.Entry(frm, width=10)
xdet_ent.insert(0, str(preset_values['xdet_ent']))
xdet_ent.grid(column = 3, row = 4)

# detail radius entry
det_radius = tk.Entry(frm, width=10)
det_radius.insert(0, str(preset_values['det_radius']))
det_radius.grid(column = 3, row = 5)

# Source-Skin Distance (SSD) entry (source-object)
ssd_ent = tk.Entry(frm, width=10)
ssd_ent.insert(0, str(preset_values['ssd_ent']))
ssd_ent.grid(column = 3, row = 6)

# Source-Detector Distance (SDD) entry
sdd_ent = tk.Entry(frm, width=10)
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
object_lst = ['tissue']                                             # object composition, although it is one material, it has been made as a list so users can further add more materials
detail_lst = ['water', 'bone', 'calcium', 'iodinated water']        # detail composition list
R_lst = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]         # scatter to primary ratio list

# tk.StringVar(parent) returns a StringVar (control variable) instance for tk.OptionMenu() variable
# set() function is used to assign previously defined preset values to the StringVar created
# tk.OptionMenu (parent, variable, list of choices), the parent is the frame, variable should be StrVariable instance and the list of choices are 

''' option menus for Phantom and Geometrical Properties '''

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

############################################### main function  #######################################################

def calculate():
       
    # making the img_32bit, contrast, field and matrix size global variables to use them in other functions.
    global img_32bit, contrast, field_of_view, matrix

    # entries for X-ray beam
    mono_energy_data = float(mono_energy_ent.get())     # monochromatic X-ray energy (keV)
    n0_data = float(n0_ent.get())                       # number of incident photons per cm²

        
    # entries for geometry and phantom properties
    sdd = float(sdd_ent.get())                  # source-detector distance
    ssd = float(ssd_ent.get())                  # source-skin (object) distance
    obj_material = obj_ent.get()                # object type
    det_material = det_ent.get()                # detail type 
    det_rad_cm = float(det_radius.get())        # detail radius in cm
    xobj = float(xobj_ent.get())                # object thickness
    xdet = float(xdet_ent.get())                # detail thickness

    # previously we have defined all_mu dictionary.
    # now we define object and detail attenuation coefficient dictionaries based on the object and detail material and energy selected by user in GUI.
    obj_dict = all_mu[obj_material][mono_energy_data]
    det_dict = all_mu[det_material][mono_energy_data]
    
    # entries used for image formation
    pixel = float(pixel_size.get())                 # piexel size
    scatter_ratio = float(scatter_ent.get())        # scatter to primary ratio
    field_of_view = float(FOV_ent.get())             # field of view
    
    # basic calculations for image formation
    matrix = int(field_of_view/pixel)               # defining number of pixels in every direction of field of view
    det_rad_pixel = det_rad_cm/pixel                # converting detail radius from cm to pixel for image formation  
    
    #################################################### Errors #############################################################

    # show error if energy is not between 10-101 keV:
    if mono_energy_data < 10 or mono_energy_data > 100:
        messagebox.showerror('Energy Error', 'Enter energy between 10 and 100 keV')

    # setting condition that object thickness must be greater than detail thickness
    if xdet >= xobj:
        messagebox.showerror('Geometry Error', 'Enter object thickness larger than detail thickness')
    
    # setting condition that SDD (source detector distance) must be greater than SSD (source skin distance)
    if ssd >= sdd:
        messagebox.showerror('Geometry Error', 'Enter source skin distance smaller than source detector distance')

    # setting condition that FOV (field of view) must be greater than detail radius to form the image.
    if det_rad_cm >= field_of_view:
        messagebox.showerror('Geometry Error', 'Enter field of view larger than detail radius')

    ############################### fluence of X-ray spectrum (photons/cm²) after attenuation ##################################

    # calculate number of photons using Beer-Lambert's law
    # obj_dict and det_dict give the linear attenuation coefficient based on the energy and material selected by user
    # np.exp : use numpy exponentional function for calculation
    # n_out = number of photons per cm² outside the detail, which degrades exponentionally (Beer-Lambert's law)
    # n_in = number of photons per cm² passing through the detail
    n_out = n0_data *np.exp(-obj_dict*xobj)                                # photons/cm² in the background
    n_in = n0_data *np.exp((-obj_dict*(xobj-xdet))-(det_dict*(xdet)))      # photons/cm² behind the detail
        
    ###################################### Rescaling using Inverse square law ############################################
                     
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

# define show_contrast function, so when the user clicks on contrast button, it recalls this function.
def show_contrast():
    # set title and messge of the message box, use f sting method to show txt and show the contrast value to 4 decimal places
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
