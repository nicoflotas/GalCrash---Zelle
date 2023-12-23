"""
CS152 - Computational Thinking: Science
Nico Flota Sanchez
08/December/2023

Project_10

Creates a galaxy merger simulation.

To run this program, simply click on the run button. Provide the parameters from the entry boxes, enter return, and once you've inputted all parameters, click the start button.
To quit, click the exit button inside the GUI.
"""

from graphicsPlus import *
from astro_objects import *

#Simulation box ranges: 
# x = (570, 1290)
# y = (5, 740) 
# or Point(570, 5), Point(1290, 740)


def main():
    """Main"""
    # create a graphics window
    win = GraphWin("Galaxy Merger", 1300, 750, False)
    win.setBackground("black")

    #create the scene
    buildSeparations(win)
    buildText(win)
    entries = buildEntryBoxes(win) #draws the entry boxes

    #get starting parameters from entry boxes
    yellowtheta = int(entries[0].get_value()) #comp
    yellowphi = int(entries[1].get_value()) #comp
    bluetheta = int(entries[2].get_value()) #gal
    bluephi = int(entries[3].get_value()) #gal
    peri = int(entries[4].get_value())
    massrat = int(entries[5].get_value())
    starnum = int(entries[6].get_value())

    #Buttons
    start = Button(win, 500, 350, 60, "Start")
    pause = Button(win, 500, 380, 60, "Pause")
    exit = Button(win, 500, 410, 60, "Exit")

    #button actions
    if start.button_action(win) == True:

        while True:
            sim = Simulate(win, yellowtheta, yellowphi, bluetheta, bluephi, peri, massrat, starnum)
            #win.getMouse()
            #sim.update_plot(0.04)
            #time.sleep(0.04) # have the animation go at the same speed
            win.update()

            """if pause.button_action(win) == True: 
                break #interrups the while loop""" #lacks functionality without an animation

            if exit.button_action(win) == True: 
                break #closes the program by ending the loop
    
    """if pause.button_action(win) == True: 
                #win.close() #exits the program"""
    
    win.close()
    
def buildSeparations(win):
    """Draws separations in the window and the explanation box."""
    separations = []

    user_box = Rectangle(Point(10, 5), Point(560, 440))
    user_box.setOutline("white")
    separations.append(user_box)

    simulation_box = Rectangle(Point(570, 5), Point(1290, 740))
    simulation_box.setOutline("white")
    separations.append(simulation_box)

    text_box = Rectangle(Point(10, 450), Point(560, 740))
    text_box.setOutline("white")
    separations.append(text_box)

    for sep in separations:
        sep.draw(win)


def buildText(win):
    """Write the instructions on the text window."""
    
    angles = Text(Point(270, 50), "Blue refers to the main galaxy, and yellow to its companion galaxy. \n Use numbers between -90 to 90 for these parameters.")
    angles.setTextColor("white")
    angles.draw(win)
    
    text = Text(Point(285, 480), "Recommended parameters:")
    text.setTextColor("white")
    text.setStyle("bold")
    text.draw(win)
    
    peritxt = Text(Point(200, 550), "Peri (distance of closest approach) = 50 kpc")
    peritxt.setTextColor("white")
    peritxt.draw(win)   

    massratxt = Text(Point(260, 600), "Mass ratio: for major mergers, set this = 1. For minor mergers, \n set this = 8. Variation exists in-between.")
    massratxt.setTextColor("white")
    massratxt.draw(win)

    numstartxt = Text(Point(260, 650), "Number of stars will be divided equally between both galaxies.")
    numstartxt.setTextColor("white")
    numstartxt.draw(win)

def buildEntryBoxes(win):
    """Gets the parameters needed to build the simulation."""

    entries = []

    #Entry box texts
    yellow_theta_text = Text(Point(100, 150), "Yellow Theta:")
    yellow_theta_text.setTextColor("white")
    yellow_theta_text.draw(win)

    yellow_phi_text = Text(Point(300, 150), "Yellow Phi:")
    yellow_phi_text.setTextColor("white")
    yellow_phi_text.draw(win)

    blue_theta_text = Text(Point(100, 250), "Blue Theta:")
    blue_theta_text.setTextColor("white")
    blue_theta_text.draw(win)

    blue_phi_text = Text(Point(300, 250), "Blue Phi:")
    blue_phi_text.setTextColor("white")
    blue_phi_text.draw(win)
    
    peri_text = Text(Point(80, 350), "Peri [kpc]:")
    peri_text.setTextColor("white")
    peri_text.draw(win)

    massrat_text = Text(Point(200, 380), "Blue Galaxy Mass (w.r.t. Yellow Galaxy Mass):")
    massrat_text.setTextColor("white")
    massrat_text.draw(win)

    numstar_text = Text(Point(100, 410), "Number of stars:")
    numstar_text.setTextColor("white")
    numstar_text.draw(win)


    #Entry boxes
    yellow_theta = EntryBox(win, 170, 150, 3)
    entries.append(yellow_theta)

    yellow_phi = EntryBox(win, 370, 150, 3)
    entries.append(yellow_phi) #recovers the text entered into the entry box

    blue_theta = EntryBox(win, 170, 250, 3)
    entries.append(blue_theta)

    blue_phi = EntryBox(win, 370, 250, 3)
    entries.append(blue_phi)

    peri = EntryBox(win, 140, 350, 3)
    entries.append(peri)

    massrat = EntryBox(win, 380, 380, 3)
    entries.append(massrat)

    numstar = EntryBox(win, 180, 410, 3)
    entries.append(numstar)


    return entries


def buildSimulation(win):
    """Creates the scene from the parameters given."""
    shapes = []

    return shapes

if __name__ == "__main__":
    main()