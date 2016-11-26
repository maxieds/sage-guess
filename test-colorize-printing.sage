from ColorizeString import *

print_color_string("Some text ... \n", bold = True, foreground = 'FGBLUE')
print_color_string("Some text ... \n", bold = True, foreground = 'FGBLUE', 
                   background = 'BGBLACK', ul = True, ol = True)

fgcolors = ['FGBLACK', 'FGRED', 'FGGREEN', 'FGYELLOW', 'FGBLUE', 
            'FGMAGENTA', 'FGCYAN', 'FGWHITE', 'FGNORMAL']
bgcolors = ['BGBLACK', 'BGRED', 'BGGREEN', 'BGYELLOW', 'BGBLUE', 
            'BGMAGENTA', 'BGCYAN', 'BGWHITE', 'BGNORMAL']
font_effects = ['BOLD', 'IT', 'UL', 'OL', 'FRAME', 'ENCIRCLE', 'INVERT', 'ST']

for fg in fgcolors: 
     for bg in bgcolors: 
          print_color_string('%s / %s | ' % (fg, bg), 
                             foreground = fg, background = bg)
          for fe in font_effects: 
               print_color_string('%s | ' % fe, opts = [fe], 
                                  foreground = fg, background = bg)
          ##
          print_color_string('\n', normal_mode = True)
     ##
##




