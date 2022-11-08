import re,sys, argparse
import matplotlib.pyplot as plt


#create a class for the string to allow the position to be navigated
#also designed for rapid analysis of trees without conversion to hierarchical
#structure
class newick_string:
    
    def __init__(self):
        self.string = ""
        self.pos = 0
        self.pos2 = 0
        self.label = []
        self.biparts = []
        self.branch_lengths = []
        self.tip_names = []
        self.tip_lengths = []
        #self.tree = Node()
    
    #designed to work with get tree interior and 
    def get_info(self,newick):
        
        name_array = []
        count = 0
        brlen = ""
        while self.pos2 < len(newick):
            #print(newick[self.pos2])
            while newick[self.pos2] == "(":
                count += 1
                self.pos2 += 1
                name = ""
            if newick[self.pos2 - 1] != ")":
                while re.match(r"^[A-Za-z0-9@\'\._-]*$", newick[self.pos2]):
                    name += newick[self.pos2]
                    self.pos2 += 1
            else:
                while re.match(r"^[A-Za-z0-9@\'\._-]*$", newick[self.pos2]):
                    self.pos2 += 1
            if newick[self.pos2] == ":":
                brlen = ""
                self.pos2 += 1
                while re.match(r"^[0-9\.:e-]*$", newick[self.pos2]):
                    brlen += newick[self.pos2]
                    self.pos2 += 1
            if newick[self.pos2] == ",":
                self.pos2 += 1
                if name != "":
                    name_array.append(name)
                    if len(self.biparts) == 0:
                        self.tip_names.append(name)
                        self.tip_lengths.append(brlen)
                name = ""
            if newick[self.pos2] == ")":
                if name != "":
                    name_array.append(name)
                    if len(self.biparts) == 0:
                        self.tip_names.append(name)
                        self.tip_lengths.append(brlen)
                    name = ""
                        
                self.pos2 += 1
                count -= 1
                #This would be the bracket that matches the top
                if count == 0:
                    label = ""
                    while re.match(r"^[A-Za-z0-9@\'\._-]*$", newick[self.pos2]):
                        label += newick[self.pos2]
                        self.pos2 += 1
                    if newick[self.pos2] == ":":
                        brlen = ""
                        self.pos2 += 1
                        while re.match(r"^[0-9\.:e-]*$", newick[self.pos2]):
                            brlen += newick[self.pos2]
                            self.pos2 += 1
                    
                    #get the other side of the bipartition
                    other_side = []
                    bip_array = []
                    if len(self.biparts) != 0:
                        #get the difference between the first which has all tips and the current
                        other_side = [x for x in self.biparts[0][0] if x not in name_array]
                    bip_array.append(name_array)
                    bip_array.append(other_side)
                    self.biparts.append(bip_array)
                    self.branch_lengths.append(brlen)
                    self.label.append(label)
                    self.pos2 = 0
                    break

    #anytime there is an open bracket grab all info
    def get_tree_interior(self):
    
        #main loop
        while self.pos < len(self.string):
            #At every bracket get the info
            if self.string[self.pos] == "(":
                self.pos2 = 0
                self.get_info(self.string[self.pos:])
            self.pos += 1
        #remove the first element which will just have everything
        #self.biparts.pop(0)
        #self.label.pop(0)
        #self.branch_lengths.pop(0)



def get_matches(tree1,tree2):

    HASH = {}
    for x in range(0,len(tree1.biparts)):
        
        for y in range(0,len(tree2.biparts)):
        
            if sorted(tree1.biparts[x][0]) == sorted(tree2.biparts[y][0]): 
            #or sorted(tree1.biparts[x][1]) == sorted(tree2.biparts[y][0]):
                
                #print str(tree1.label[x]) + "\t" + str(tree2.label[y])
                HASH[tree1.label[x]] = tree2.label[y]
                #print str(tree1.biparts[x]) + "\t" + str(tree1.label[x]) + "\t" + str(tree2.label[x])
                #break
    return HASH         
            
def get_csv_details(analysis_file,breaks):
    
    InfoHash = {}
    analysis_file.readline() #hacky way to remove header
    
    for line in analysis_file:
        
        line_array = []
        line = line.strip("\n\r")
        line_array = line.split(",") #pos 1 matched node, pos 4 number of conflicts
        if line_array[0] not in InfoHash:
            InfoHash[line_array[0]] = []
        InfoHash[line_array[0]].append(line_array[4])
    return InfoHash

        

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="A program to draw a pie chart from the output of conflict_detector.py.")
    parser.add_argument("-f", "--folder", type=str, help="The location of all the input files - can be a full path or just from current working directory.")
    parser.add_argument("-p", "--input_prefix", type=str, help="The prefix of input files (X, where X is X_labels.tre, X_concord.tre, etc).")
    parser.add_argument("-t", "--totals_tree", type=str, help="The *_total_analyzed.tre file that corresponds.")
    parser.add_argument("-a", "--alternatives_reported", type=int, help="The number of alternative conflicts to be reported, anything up to 15.")
    parser.add_argument("-o", "--output_type", type=str, default="n", help="Should the pie chart be labelled (l), emphasized (e), both (le) or neither (default)?")
    parser.add_argument("-at", "--alt_total", type=int, help="This program uses the total number of times a node was analysed as the total to draw pie charts, but if you are using ortholog trees you may want to use the total number of genes that were included in the analysis. If that is the case, you can input the total number of genes here.")

    if len(sys.argv) == 1:
        parser.print_usage()
        parser.exit()

    args = parser.parse_args()

    if args.folder[-1] == "/":
        folder = args.folder
    else:
        folder = args.folder + "/"
    prefix = args.input_prefix
    totals = args.totals_tree
    breaks = args.alternatives_reported
    out_type = args.output_type

    #Steps
    #1) Read in the tree with annotations
    name1 = str(folder) + str(prefix) + "_labels.tre"
    tree_file = open(name1, "r")
    tree1 = newick_string()
    tree1.string = tree_file.readline().strip("\n")
    tree1.get_tree_interior()
    
    #2) Read in the total analyzed tree
    name2 = str(totals)
    tree_file2 = open(name2, "r")
    tree2 = newick_string()
    tree2.string = tree_file2.readline().strip("\n")
    tree2.get_tree_interior()

    #3) Read in tree with concordant
    name3 = str(folder) + str(prefix) + "_concord.tre"
    tree_file3 = open(name3, "r")
    tree3 = newick_string()
    tree3.string = tree_file3.readline().strip("\n")
    tree3.get_tree_interior()
        
    #4) Read in tree with conflict
    name4 = str(folder) + str(prefix) + "_conflict.tre"
    tree_file4 = open(name4, "r")
    tree4 = newick_string()
    tree4.string = tree_file4.readline().strip("\n")
    tree4.get_tree_interior()
    
    #5) Read in the *_analysis.csv file
    name5 = str(folder) + str(prefix) + "_analysis.csv"
    analysis_file = open(name5, "r")
    
    #Yay next part
    #Steps
    #1 Get the matches between the annotation and the other trees
    #2 Read in the analysis file and get the number of breaks
    #3 Turn into pie charts
    
    #1) Get the matches between the annotation and the other trees
    TotalAnalyzedHash = get_matches(tree1,tree2)
    TotalConcordantHash = get_matches(tree1,tree3)
    TotalConflictHash = get_matches(tree1,tree4)
    
    #2 Read in the analysis file and get the number of breaks
    HashOfAltCons = {}
    HashOfAltCons = get_csv_details(analysis_file,breaks)
    
    #3 Turn into pie charts
    #Ugh...
    #Piecharts are the concordant value, the alternative conflicts based on breaks specified
    #the conflict value minus the sum of the alternative conflicts, the difference between everything
    #combined and the total analyzed is the uninformative
    
    
    #alternative colors
    colors = ["yellow", "orange", "green", "purple", "sienna", "magenta", "pink", "burlywood", "darkseagreen", "teal", "turquoise", "honeydew", "palegoldenrod", "ivory", "rosybrown"]
    hatches = ['//', '\\\\', 'xx', '///', '\\\\\\', 'xxx']


    for number in range(0,len(TotalAnalyzedHash)):
        
        concordant_value = TotalConcordantHash[str(number)]
        conflict_value = int(TotalConflictHash[str(number)])
        if args.alt_total:
            total_uninformative = int(args.alt_total)
        else:
            total_uninformative = int(TotalAnalyzedHash[str(number)]) #will have necessary values subtracted
        
        #check that there actually is an alternative conflict
        if str(number) in HashOfAltCons:
            
            labels = ["Concordant"]
            #check that there are enough breaks
            
            if len(HashOfAltCons[str(number)]) < int(breaks):
                
                breaks = len(HashOfAltCons[str(number)])
            
            
            if breaks != 0:
                
                out_colors = ["blue"] #smaller array for what will be printed
                out_hatches = ["/"]
                out_alts = [] #will have numbers associated with alternative conflicts
                total_con = int(conflict_value) #way to maintain total conflicts
                
                #get numbers associated with conflicts
                # get suitable amount of alt conflict colours

                for numb in range(0,int(breaks)):
                    out_colors.append(colors[numb])
                    out_hatches.append(hatches[numb])
                    out_alts.append(int(HashOfAltCons[str(number)][numb]))
                    conflict_value -= int(HashOfAltCons[str(number)][numb])
                    ConName = "Conflict_" + str(numb)
                    labels.append(ConName)
            
            labels.append("OtherConflict")
            labels.append("Uninformative")
            total_uninformative = total_uninformative - int(total_con) - int(concordant_value) 
            
            #make it match the colors array
            size_array = []
            size_array.append(int(concordant_value))    
            size_array += out_alts
            size_array.append(conflict_value)
            size_array.append(total_uninformative)
            
            #put in the red and grey
            out_colors.append("red")
            out_colors.append("silver")
            out_hatches.append("\\")
            out_hatches.append("x")
            
            # print size_array
            # print out_colors
            # print out_hatches

            fig1, ax1 = plt.subplots(1, 1, squeeze = False)
            if out_type == "l":
                piechart = ax1[0,0].pie(size_array, colors = out_colors, labels = labels)
            
            elif out_type == "le":
                explode = []
                explode.append(0.2)
                for x in range(0,(len(out_colors)-1)):
                    explode.append(0.1)
                piechart = ax1[0,0].pie(size_array, colors = out_colors, labels = labels, explode = explode)

            elif out_type == "e":
                explode = []
                explode.append(0.2)
                for x in range(0,(len(out_colors)-1)):
                    explode.append(0.1)
                piechart = ax1[0,0].pie(size_array, colors = out_colors, explode = explode)
            else:   
                piechart = ax1[0,0].pie(size_array, colors = out_colors)

            # crosshatching on the pie charts 

            ax1[0,0].axis('equal')
            
            for i in range(len(piechart[0])):
                piechart[0][i].set_hatch(out_hatches[(i)])
            


            #plt.show()
            fig_name = "Node_" + str(number) + ".svg"
            plt.savefig(fig_name)
            print("\n" + "Piechart For Node " + str(number) + "\nRelationship\tColor")
            for x in range(0,len(out_colors)):
                print(labels[x] + "\t" + out_colors[x])
            
        #So all is concordant or uninformative
        else:
            print("No Conflict for node: " + str(number))

    

    #labels = ["Frogs", "Hogs", "Dogs", "Logs"]
    #sizes = [15, 30, 45, 100]
    #explode = (0, 0.1, 0, 0)  # only "explode" the 2nd slice (i.e. 'Hogs')

    #fig1, ax1 = plt.subplots()
    #ax1.pie(sizes, explode=explode, labels=labels, autopct='%1.1f%%',
    #    shadow=True, startangle=90)
    #ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
    
    #test = "1" + ".svg"
    #plt.savefig(test)
    #plt.show()
