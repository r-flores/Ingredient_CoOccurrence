#Author: Ricky Flores
#Orginal Work: https://github.com/kmcooper/Ingredient_CoOccurrence
#Created: 1/3/20
#Last editied: 1/6/20
#input: ingredients.tab
#output: ingredient frequency from a hash
#notes: this is a python adaptation to https://github.com/kmcooper/Ingredient_CoOccurrence/blob/master/code/ingredient_network.pl

import getopt
import sys
import os
import re

# Function Name: Usage
# Parameters: n/a
# Return: Warning message about parameters
def usage():
    print(("\nUsage: ingredient_has.py [-i ingredients.tab] [-o FILE]")
    + "\n \t-i: ingredients.tab = The file generated by the Rscript"
    + "\n \t-o: outputfile.output = name fo the file the network will be written to")

infile = ""
outfile = ""

# Create a hash for all ingredients
edge_hash = {}

try:
    opts, args = getopt.getopt(sys.argv[1:], 'i:o:')
except getopt.GetoptError as err:
    sys.stdout = sys.stderr
    print(str(err))
    usage()
    sys.exit(2)

for (opt, args) in opts:
    if (opt == '-i'):
        print("Received the following file as input: " + args)
        infile = args
    if (opt == '-o'):
        print("Received the following file as network output: " + args)
        outfile = args

# Desc.: routine for adding items to the hash
# Input: ingredient array (each ingredient is called an item)
# return: None
def make_food_subnetwork(i):
    foods = i
    foods.sort()
    arrsize = len(foods)
    for i in range(len(foods)):
        arrsize -= 1
        tempi = i
        for j in range(arrsize):
            node1 = foods[i]
            node2 = foods[tempi + 1]
            node1 = re.sub(r"\s", "_", node1)
            node2 = re.sub(r"\s", "_", node2)
            tempi += 1
            key = node1 + " " + node2
            if (key in edge_hash):
                edge_hash[key] += 1
            else:
                edge_hash[key] = 1
        

# number of foods
dataLength= open(infile, encoding = "utf8")
foodCount = len(dataLength.readlines())
dataLength.close

# open the file
data = open(infile, encoding = "utf8")
line = data.readline()

# An array for bulding co-occurence
ing_arr = []

# count of lines
countline = 0

# count of foods in a single line
single_ingre_foods = 0

# Iterate through each line
# Each line is one food with multiple ingredients
while (line):

    # inc the number of lines for each iteration
    countline += 1

    # save the current current line
    currentLine = line

    # print progress
    #print ("Processing food " + str(countline) + " of " + str(foodCount))

    #Variables for barcode and ingredient
    #barcode should always be positive
    barcode = -10000
    ingredients = ""

    #split the line on tab between barcode and ingredients
    tabs = currentLine.split("\t")

    #split the line on tab between barcode and ingridents
    tabs = currentLine.split("\t")
    if (len(tabs) >= 1):
        barcode = tabs[0]
    else:
        #call out lines with no barcode
        print("no barcode found on line " + str(countline))

    if (len(tabs) > 1):
        ingredients = tabs[1]

        # remove newline, replace "and", and "contains" with commas
        ingredients = ingredients.replace("\n", "")
        ingredients = ingredients.replace("(", ",")
        ingredients = ingredients.replace(")", ",")
        ingredients = ingredients.replace(" and ", ",")
        ingredients = ingredients.replace(" contains ", ",")

        #seperate ingredients by commas
        #check if current ingredient has commas
        if (re.search(",", ingredients)):
            
            #split ingredients by comma to sperate each item
            items = ingredients.split(",")
            
            #iterate through each item
            for item in items:
                curr_item = item
                curr_item = re.sub(r'^\s+|\s+$', '', curr_item)
                curr_item = re.sub(r'\(.*\)', '', curr_item)
                curr_item = re.sub(r'\(.*', '', curr_item)
                curr_item = re.sub(r'.*\)', '', curr_item)

                # if current item is empty
                if (re.search(r"^\s+$", curr_item)):
                    #Do nothing
                    pass

                # if current item is only one word
                elif (re.search(r"^\s*(\w+)\s*$", curr_item)):
                    item = curr_item
                    ing_arr.append(item)

                # if current item contains and, with, or, and it should be a comma
                elif (re.search(r".*( and | with | or )+.*", curr_item)):
                    splitonPattern = re.compile(r".*( and | with | or )+.*")
                    spliton = splitonPattern.search(curr_item).groups()[0]
                    items = curr_item.split(spliton)
                    for item in items:
                        citem = item
                        citem = re.sub(r'^\s+', '', citem)
                        citem = re.sub(r'\s+$', '', citem)
                        if (bool(re.search(r"^\s+$", citem)) == False):
                            ing_arr.append(item)
                        else:
                            print("Error: unable to process this item: " + citem)

                # if items are multiple words    
                elif (re.search(r"\s*([\S+\s+]+)\s*", curr_item)):
                    item = curr_item
                    ing_arr.append(item)

                # If the current item is an empty string
                elif (re.search(r"^(?![\s\S])", curr_item)):
                    # The item is empty so move on
                    pass
                
                ###############################################
                ### Lines omitted from perl script version  ###
                ###############################################  
                else:
                    print("Unable to process item: " + curr_item)

        # no commas
        else:

            # remove empty parentheses
            ingredients = re.sub(r'\(.*\)', '', ingredients)
            ingredients = re.sub(r'\(', '', ingredients)
            ingredients = re.sub(r'\)', '', ingredients)

            # if ingredients are empty space
            if (re.search(r"^\s+$", ingredients)):
                # do nothing
                pass

            # if ingredients are just one word
            elif (re.search(r"^\s*(\w+)\s*$", ingredients)):
                item = ingredients
                ing_arr.append(item)

            # if ingredients contain and, with, or or
            elif (re.search(r".*( and | with | or )+.*", ingredients)):
                #splitonPattern = re.compile(r".*( and | with | or )+.*")
                #spliton = splitonPattern.search(curr_item).groups()[0]
                spliton = re.findall(r".*( and | with | or )+.*", ingredients)[0]
                items = curr_item.split(spliton)
                for item in items:
                    citem = item
                    ing_arr.append(citem)

            # if ingredients contain multiple words
            elif (re.search(r"\s*([\S+\s+]+)\s*", ingredients)):
                item = ingredients
                ing_arr.append(item)

            # if the current item is an empty string
            elif (re.search(r"^(?![\s\S])", ingredients)):
                #empty item no space simply move on
                pass

            else:
                print ("Error: unable to porcess this item: " + ingredients)
    else:
        print("No ingredients found on line " + str(countline))

    # check if the array has atleast 2 entries
    if(len(ing_arr) > 1):
        make_food_subnetwork(ing_arr)
    else:
        if ing_arr:
            print (ing_arr[0] + " Not added as it was the only ingredient listed for food " + barcode)
            single_ingre_foods += 1

    ing_arr = []
    #move to the next line (food item)
    line = data.readline()
data.close

print ("A total of " + str(single_ingre_foods) + " single ingredient foods were found and are not represented in the network\n")

out = open(outfile, "w", encoding='utf-8')
for key in edge_hash:
    toWrite = (key + " " + str(edge_hash[key]) + '\n')
    out.write(toWrite)
out.close
