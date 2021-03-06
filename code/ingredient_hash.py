#Author: Ricky Flores
#Original Work: https://github.com/kmcooper/Ingredient_CoOccurrence
#Created: 11/25/19
#Last editied: 1/6/20
#input: ingredients.tab
#output: ingredient frequency from a hash
#notes: this is a python adaptation to https://github.com/kmcooper/Ingredient_CoOccurrence/blob/master/code/ingredient_hash.pl

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
hashTable = {}

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

#Function Name: add_item
#parameters: ingredient item
#Return: n/a
def add_item(i):
    i = i.replace('\n', '')
    i = re.sub(r'^\s+|\s+$', '', i)
    if (i in hashTable):
        hashTable[i] += 1
    else:
        hashTable[i] = 1

# open the file
data = open(infile, encoding = "utf8")
line = data.readline()

# count of lines
countline = 0

# iterate through each line
# each line is one food with multiple ingreidents
while (line):

    # inc the number of lines for each iteration
    countline += 1
    #save the current current line
    currentLine = line

    #Variables for barcode and ingredient
    #barcode should always be positive
    barcode = -10000
    ingredients = ""

    #split the line on tab between barcode and ingridents
    tabs = currentLine.split("\t")
    if (len(tabs) >= 1):
        barcode = tabs[0]
    else:
        #call out lines with no barcode
        print("no barcode found on line " + str(countline))

    # putting ingredients into the hash
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
                    add_item(item)

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
                            add_item(citem)
                        else:
                            print("Error: unable to process this item: " + citem)

                # if items are multiple words    
                elif (re.search(r"\s*([\S+\s+]+)\s*", curr_item)):
                    item = curr_item
                    add_item(item)

                # If the current item is an empty string
                elif (re.search(r"^(?![\s\S])", curr_item)):
                    # The item is empty so move on
                    pass
                
                #########################################################
                ### Lines 170 - 201 omitted from perl script version  ###
                #########################################################  
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
                add_item(item)

            # if ingredients contain and, with, or or
            elif (re.search(r".*( and | with | or )+.*", ingredients)):
                #splitonPattern = re.compile(r".*( and | with | or )+.*")
                #spliton = splitonPattern.search(curr_item).groups()[0]
                spliton = re.findall(r".*( and | with | or )+.*", ingredients)[0]
                items = curr_item.split(spliton)
                for item in items:
                    citem = item
                    add_item(citem)

            # if ingredients contain multiple words
            elif (re.search(r"\s*([\S+\s+]+)\s*", ingredients)):
                item = ingredients
                add_item(item)

            # if the current item is an empty string
            elif (re.search(r"^(?![\s\S])", ingredients)):

                #empty item no space simply move on
                pass

            else:
                print ("Error: unable to porcess this item: " + ingredients)
    else:
        print("No ingredients found on line " + str(countline))
    #move to the next line (food item)
    line = data.readline()
data.close

out = open(outfile, "w", encoding='utf-8')
SortedHashTable = sorted(hashTable, key=lambda i: int(hashTable[i]))

print ("Writing frequency to hash: " + outfile)
ingredientCount = 0
for key in SortedHashTable:
    toWrite = (key + " --> " + str(hashTable[key]) + "\n")
    ingredientCount += 1
    out.write(toWrite)

print(str(ingredientCount) + " ingredients found total")
out.close()
