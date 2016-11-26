#### ColorizeString.py : 
#### Helper functions to colorize terminal-printed strings by ANSI 
#### escape code colors and related formatting
#### Author:  Maxie D. Schmidt
#### Created: 2016.11.23

import sys

ansi_format_dict_init = {
     'FGBLACK'   : 30,  # Foreground colors
     'FGRED'     : 31, 
     'FGGREEN'   : 32, 
     'FGYELLOW'  : 33, 
     'FGBLUE'    : 34, 
     'FGMAGENTA' : 35, 
     'FGCYAN'    : 36, 
     'FGWHITE'   : 37, 
     'FGNORMAL'  : 39, 
     'BGBLACK'   : 40,  # Background colors
     'BGRED'     : 41, 
     'BGGREEN'   : 42, 
     'BGYELLOW'  : 43, 
     'BGBLUE'    : 44, 
     'BGMAGENTA' : 45, 
     'BGCYAN'    : 46, 
     'BGWHITE'   : 47, 
     'BGNORMAL'  : 49, 
     'BOLD'      : 1,   # Bold text 
     'IT'        : 3,   # Italicized text 
     'UL'        : 4,   # Underlined text
     'OL'        : 53,  # Overlined text 
     'FRAME'     : 51,  # Framed text
     'ENCIRCLE'  : 52,  # Encircled text
     'INVERT'    : 7,   # BG/FG Inverted text
     'ST'        : 9,   # Strikethrough text 
     'NORMAL'    : 0,   # Set back to normal mode 
     ''          : 0, 
}; 
ansi_format_dict = dict(ansi_format_dict_init) 

def get_escape_code(options): 
     format_str = '\033[0'
     for opt in options: 
          if opt not in ansi_format_dict: 
               raise NameError('Colorize option \'' + opt + '\' not supported.')
          ##
          format_str += ';%d' % ansi_format_dict[opt]
     ## 
     format_str += 'm'
     return format_str
## def 

def colorize_string(s, opts = [], print_string = False, 
                    foreground = None, background = None, 
                    ul = None, ol = None, framed = None, encircled = None, 
                    bold = None, invert = None, strikethrough = None, 
                    normal_mode = None): 

     if foreground != None:
          opts += [foreground]
     if background != None: 
          opts += [background] 
     ## 

     fnparams = [
          (ul, 'UL'), (ol, 'OL'), (framed, 'FRAME'), (encircled, 'ENCIRCLE'), 
          (bold, 'BOLD'), (invert, 'INVERT'), (strikethrough, 'ST'), 
          (normal_mode, 'NORMAL'), 
     ]
     for (param, pstr) in fnparams: 
          if param != None: 
               opts += [pstr]
          ##
     ##
     escape_code = get_escape_code(opts)
     normal_mode = get_escape_code(['NORMAL'])
     rstr = escape_code + s + normal_mode
     if print_string: 
           sys.stdout.write(rstr)
           sys.stdout.flush()
     ## 
     return rstr

## def 

def print_color_string(s, opts = [], foreground = None, background = None, 
                       ul = None, ol = None, framed = None, encircled = None, 
                       bold = None, invert = None, strikethrough = None, 
                       normal_mode = None): 

     colorize_string(s, opts = opts, foreground = foreground, 
                     background = background, 
                     ul = ul, ol = ol, framed = framed, encircled = encircled, 
                     bold = bold, invert = invert, 
                     strikethrough = strikethrough, normal_mode = normal_mode, 
                     print_string = True)

## def 

