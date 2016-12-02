
r"""
ColorizeString.py 

String processing and helper functions that add ANSI escape code 
colors and formatting to the string for pretty printing 
in a terminal. 

AUTHORS: 
- Maxie D. Schmidt (Created: 2016.11.23) 

""" 

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
     r""" 
     Returns the ANSI escape code corresponding to input formatting options. 

     INPUT: 
     - ``options`` -- a list of constants specifying formatting and 
                      color options. Each option in the list should be one of the 
                      following entries in the ansi_format_dict defined in the 
                      source above: 
                      'FGBLACK', 'FGRED', 'FGGREEN', 'FGYELLOW', 
                      'FGBLUE', 'FGMAGENTA', 'FGCYAN', 'FGWHITE', 'FGNORMAL', 
                      'BGBLACK', 'BGRED', 'BGGREEN', 'BGYELLOW', 
                      'BGBLUE', 'BGMAGENTA', 'BGCYAN', 'BGWHITE', 'BGNORMAL', 
                      'BOLD', 'IT', 'UL', 'OL', 'FRAME', 'ENCIRCLE', 
                      'INVERT', 'ST', 'NORMAL'. 

     EXAMPLES: 
     Create a bold underlined line of text with foreground color black and 
     background color yellow: 
     sage: from ColorizeString import *
     sage: s = get_escape_code(['FGBLACK', 'BGYELLOW', 'BOLD', 'UL'])
     """
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
     r""" 
     Returns a colorized and formatted string according to the boolean-valued 
     options passed to the function. 

     INPUT: 
     - ``s`` -- The string text to be formatted 
     - ``opts`` -- An optional list of additional formatting options 
                   (see the allowed options in the documentation for 
                   get_escape_code)
     - ``print_string`` -- Specifies whether the formatted string should be 
                           printed to stdout before it is returned 
     - ``foreground`` -- The foreground text color. If set, it should be one of: 
                         'FGBLACK', 'FGRED', 'FGGREEN', 'FGYELLOW', 
                         'FGBLUE', 'FGMAGENTA', 'FGCYAN', 'FGWHITE', 'FGNORMAL' 
     - ``background`` -- The background text color. If set, it should be one of: 
                         'BGBLACK', 'BGRED', 'BGGREEN', 'BGYELLOW', 
                         'BGBLUE', 'BGMAGENTA', 'BGCYAN', 'BGWHITE', 'BGNORMAL' 
     - ``ul`` -- Specifies whether the string is underlined 
     - ``ol`` -- Specifies an overline for the string
     - ``framed`` -- Whether the string should be framed (generally unsupported)
     - ``encircled`` -- Whether the string shoulde be encircled (unsupported)
     - ``bold`` -- Whether the string is bold-formatted 
     - ``invert`` -- Whether the foreground and background colors should be 
                     inverted. 
     - ``strikethrough`` -- Whether the text is struckthrough
     - ``normal_mode`` -- Whether the string is returned in normal mode 
                          with no formatting

     EXAMPLES: 
     See the documentation for print_color_string below. 
     """
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
     r""" 
     Prints a colorized and formatted string according to the boolean-valued 
     options passed to the function. 

     INPUT: 
     - ``s`` -- The string text to be formatted 
     - ``opts`` -- An optional list of additional formatting options 
                   (see the allowed options in the documentation for 
                   get_escape_code)
     - ``foreground`` -- The foreground text color. If set, it should be one of: 
                         'FGBLACK', 'FGRED', 'FGGREEN', 'FGYELLOW', 
                         'FGBLUE', 'FGMAGENTA', 'FGCYAN', 'FGWHITE', 'FGNORMAL' 
     - ``background`` -- The background text color. If set, it should be one of: 
                         'BGBLACK', 'BGRED', 'BGGREEN', 'BGYELLOW', 
                         'BGBLUE', 'BGMAGENTA', 'BGCYAN', 'BGWHITE', 'BGNORMAL' 
     - ``ul`` -- Specifies whether the string is underlined 
     - ``ol`` -- Specifies an overline for the string
     - ``framed`` -- Whether the string should be framed (generally unsupported)
     - ``encircled`` -- Whether the string shoulde be encircled (unsupported)
     - ``bold`` -- Whether the string is bold-formatted 
     - ``invert`` -- Whether the foreground and background colors should be 
                     inverted. 
     - ``strikethrough`` -- Whether the text is struckthrough
     - ``normal_mode`` -- Whether the string is returned in normal mode 
                          with no formatting

     EXAMPLES: 
     
     To print all possible combinations of the formatting / colorization 
     options to the terminal, run the following loop: 
     
     sage: from ColorizeString import *
     sage: fgcolors = ['FGBLACK', 'FGRED', 'FGGREEN', 'FGYELLOW', 'FGBLUE', 
     ....: 'FGMAGENTA', 'FGCYAN', 'FGWHITE', 'FGNORMAL']
     sage: bgcolors = ['BGBLACK', 'BGRED', 'BGGREEN', 'BGYELLOW', 'BGBLUE', 
     ....: 'BGMAGENTA', 'BGCYAN', 'BGWHITE', 'BGNORMAL']
     sage: font_effects = ['BOLD', 'IT', 'UL', 'OL', 'FRAME', 'ENCIRCLE', 
     ....: 'INVERT', 'ST']
     sage: for fg in fgcolors: 
     ....:      for bg in bgcolors: 
     ....:           print_color_string('%s / %s | ' % (fg, bg), 
     ....:                              foreground = fg, background = bg)
     ....:           for fe in font_effects: 
     ....:                print_color_string('%s | ' % fe, opts = [fe], 
     ....:                                   foreground = fg, background = bg)
     ....:           print_color_string('\n', normal_mode = True)
     """
     colorize_string(s, opts = opts, foreground = foreground, 
                     background = background, 
                     ul = ul, ol = ol, framed = framed, encircled = encircled, 
                     bold = bold, invert = invert, 
                     strikethrough = strikethrough, normal_mode = normal_mode, 
                     print_string = True)
## def 

