##This contains a set of functions that returns the a style sheet 
##for each Q Widget of interest
##

#Main Window
def QMain():
    layout = """QWidget{
                       background-color: #a8a8a8;
                       color: #000099;
                       font: bold 100%;}
                QWidget:drop-down {
                       subcontrol-origin: padding;
                       subcontrol-position: top right;
                       width: 15px;
                       border-left-width: 1px;
                       border-left-color: darkgray;
                       border-left-style: solid; /* just a single line */
                       border-top-right-radius: 3px; /* same radius as the QComboBox */
                       border-bottom-right-radius: 3px;
}

             """
    return layout

#Menu Bar
def QMenuBar():
    layout = """QMenuBar {
                     background-color: #d3d3d3}
                QMenuBar:item {
                     spacing: 5px; /* spacing between menu bar items */
                     padding: 2px 8px;
                     background: #bdbdbd;
                     border-radius: 4px;}
                QMenuBar:item:selected { /* when selected using mouse or keyboard */
                     background: #939393;}
                QMenuBar:item:pressed {
                     background: #545454;}
                """
    return layout

#Menu
def QMenu():
    layout = """QMenu {
                     background-color: #d3d3d3}
                QMenu:item {
                     spacing: 5px; /* spacing between menu bar items */
                     padding: 2px 12px;
                     background: transparent;
                     border-radius: 4px;}
                QMenu:item:selected { /* when selected using mouse or keyboard */
                     background: #ffffff;}
                QMenu:item:pressed {
                     background: #d3d3d3;}
                """
    return layout

#Menu
def DisQMenu():
    layout = """QMenu {
                     background-color: #d3d3d3;
                     color: #6495ed}
                QMenu:item {
                     spacing: 5px; /* spacing between menu bar items */
                     padding: 2px 12px;
                     background: transparent;
                     border-radius: 4px;}
                QMenu:item:selected { /* when selected using mouse or keyboard */
                     background: #ffffff;}
                QMenu:item:pressed {
                     background: #d3d3d3;}
                """
    return layout

#LineEdit
def QLineEdit():
    layout = """QLineEdit {
                     background-color:#ffffff;
                     color: #000099;}
             """
    return layout

#Push Buttons
def QPushButton():
    layout = """QPushButton{
                            background-color: #d3d3d3;
                            border-style: outset;
                            border-width: 2px;
                            border-radius: 10px;
                            border-color: black;
                            color: #000099;
                            font: 18px;
                            min-width: 4em;
                            max-width: 8em;
                            padding: 6px 6px 6px 6px;}
                QPushButton:hover {
                            background-color: #bdbdbd;}      
                QPushButton:pressed {
                            background-color: #545454;
                            border-style: inset;}
             """
    return layout

def QPushButton2():
    layout = """QPushButton{
                            background-color: #d3d3d3;
                            border-style: outset;
                            border-width: 2px;
                            border-radius: 10px;
                            border-color: black;
                            color: #000099;
                            font: bold 12px;
                            min-width: 3em;
                            max-width: 5em;
                            padding: 6px 6px 6px 6px;}
                QPushButton:hover {
                            background-color: #939393;}
                QPushButton:pressed {
                            background-color: #545454;
                            border-style: inset;}
             """
    return layout

def QPushButton3():
    layout = """QPushButton{
                            background-color: #d3d3d3;
                            border-style: outset;
                            border-width: 2px;
                            border-radius: 10px;
                            border-color: black;
                            color: #000099;
                            font: bold 12px;
                            padding: 6px 6px 6px 6px;}
                QPushButton:hover {
                            background-color: #939393;}
                QPushButton:pressed {
                            background-color: #545454;
                            border-style: inset;}
             """
    return layout



#Combo box
def QComboBox():
    layout = """QComboBox {
                           border: 1px solid gray;
                           border-radius: 3px;
                           background-color: #d3d3d3;
                           color: #000099;
                           padding: 0px 0px 0px 0px;
                           min-width: 3.5em;
                           combobox-popup: 0;}
                QComboBox:editable {
                           background: black;}
                QComboBox:hover {
                           background-color: #bdbdbd;}
                QComboBox:on { /* shift the text when the popup opens */
                           padding-top: 3px;
                           padding-left: 4px;}
                QComboBox:drop-down {
                           subcontrol-origin: padding;
                           subcontrol-position: top right;
                           width: 1.5em;
                           max-height: 8em;
                           border-left-width: 1px;
                           border-left-color: darkgray;
                           border-left-style: solid; /* just a single line */
                           border-top-right-radius: 3px; /* same radius as the QComboBox */
                           border-bottom-right-radius: 3px;}    
                QComboBox:down-arrow {
                           image: url(./icons/down.png);
                           border: 0px;
                           height: 1.5em;
                           width: 0.75em;
                           color: #000000}
             """
    return layout

#Group box
def QGroupBox():
    layout = """QGroupBox{
                          border: 1px solid #000099; 
                          padding: 0px 0px 0px 0px;
                          font-family:arial;
                          font-size:18px}                                   
                QGroupBox:title{
                          subcontrol-origin: margin;
                          subcontrol-position: top center;}
             """
    return layout

#Group box - plot gbox
def QGroupBox2():
    layout = """QGroupBox{
                          border: 1px solid #000099;
                          margin-top: 1em;
                          padding: 0px 0px 0px 0px;
                          font-family:arial;
                          font-size:18px}
                QGroupBox:title{
                          subcontrol-origin: margin;
                          subcontrol-position: top center;}
             """
    return layout

#Check Box
def QCheckBox():
    layout = """QCheckBox{
                          spacing: 5px;}                                   
             """
    return layout

#Qscroll Area
def QScrollArea():
    layout = """QScrollArea{
                          border: 1px solid #000099;
                          padding: 0px 0px 0px 0px;
                          margin: 0px 0px 0px 0px;}
             """
    return layout






