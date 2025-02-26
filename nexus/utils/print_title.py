"""
Module: print_title
-------------------

This module provides a function to print the title and version of the package.

Functions:
----------
    - print_title: Prints the title and the version of the package.
"""

def print_title(__version__) -> None:
    """
    Prints the title and the version of the package.

    Parameters:
    -----------
        __version__ (str): The version of the package.

    Returns:
    --------
        None
    """
    title = r"""
                                                                                                                                    
                :-:      --:   -=+++=-  -:     ==: ::       :-:    :+*##*=:               -===-:      -=+=:    ::---::::::          
               *@@%#-   +@@# -%@@@@@@@=*@%*: :#@@+=@@+     -@@%:  =@@@@@@@+             +%@@@@@%+   -#@@@@@*  :%@@@@@@@@@@*         
              -@@@@@@+  #@@%:%@@%+==== -%@@@#%@@#:+@@#     =@@@-  %@@%--=-             *@@@#+*%@%: =@@@#*@@@*  +##%@@@@%%#+         
              =@@@@@@@= #@@#-@@@%+=-    :*@@@@@+  +@@%:    +@@@-  +@@@@%#*:   ::::::  =@@@+   :-: :@@@#  *@@@:    :%@@%:            
              =@@@=#@@@#@@@+=@@@@@@*      #@@@#   =@@@=    *@@%:   -+#%@@@%- *%%@@@@+ +@@@:       *@@@%##%@@@=     #@@#             
              =@@% :%@@@@@@==@@@*--     :*@@@@@#: :%@@@*==*@@@+ -##+  :#@@@+ +###%%#= =@@@*: :-- :%@@@@%%@@@@*     #@@%             
              =@@%  -%@@@@# :@@@#*###*:-%@@@#%@@@- -%@@@@@@@@*  *@@@#*#@@@@-           +@@@%%%@@=-@@@#:::-@@@#     *@@%:            
              -%%+   :+##+:  =%@@@@@@%:=@@%- :#%%-  :+#%%%%*-   :%@@@@@@@#=             -*%@@@%*::%@#:    +@@*     -%%+             
                :              :-----:  :-     :       :::        -=+++=:                  :-::   :-:      ::                       
                                                                                                                                    
    """
    old_title = r'''


`7MN.   `7MF'`7MM"""YMM  `YMM'   `MP'`7MMF'   `7MF'.M"""bgd
  MMN.    M    MM    `7    VMb.  ,P    MM       M ,MI    "Y
  M YMb   M    MM   d       `MM.M'     MM       M `MMb.
  M  `MN. M    MMmmMM         MMb      MM       M   `YMMNq.
  M   `MM.M    MM   Y  ,    ,M'`Mb.    MM       M .     `MM
  M     YMM    MM     ,M   ,P   `MM.   YM.     ,M Mb     dM
.JML.    YM  .JMMmmmmMMM .MM:.  .:MMa.  `bmmmmd"' P"Ybmmd"


    '''
    print(title)
    print(f"__version__ \u279c\t {__version__}\n")
    return
