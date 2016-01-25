    

class htmlTable(object): 
    """General Table formatter. It can be used to create a grid structure
    using python syntaxis, and convert it into HTML format.
    by jhcepas@gmail.com
    modified by l.p.pryszcz@gmail.com
    ##
    # Example: 

    T = htmlTable()
    T.add_cell(0,"Phylome")
    T.add_cell(0,"Description")
    T.add_cell(1,"Human Phylome")
    T.add_cell(1,"2007 Genome Biology paper.")

    # Set CSS style used for different rows in table

    T.table_style = "gtable"
    T.header = True 

    print T.asHTML()
    """
    def __init__(self):
        self.rows = []
        self.header = True
        self.table_style = ""
        self.header_style = ""
        self.content_style = ""
        self.td_flags = ""
        self.row_classes = {}
        self.cell_classes= {}

    def add_cell(self, row_index, content, classname=None):
        """ Add the given object (any text or HTML code) as the content of
        a cell in the row specified by the "row_index" argument. Cells
        are arranged in a stack-like way. """
        # Creates the grid until fitting the given index
        while len(self.rows) <= row_index:
            self.rows.append([])
        self.rows[row_index].append(str(content))
        if classname: 
            self.cell_classes[row_index][len(self.rows[row_index]-1)] = classname
            
    def remove_column(self, column_index):
        """Remove column from table"""
        for row_index, row in enumerate(self.rows):
            if len(row)<=column_index:
                continue
            #remove column content
            row.pop(column_index)
            self.rows[row_index] = row
            #drop formatting classes
            #if self.cell_classes[row_index]:

    def asHTML(self):
        """Returns the HTML representation of the table object."""
        html = ['<table align=center class="%s">' % self.table_style]
        for r in self.rows:
            html.append(' <tr>\n  ')
            for cell in r:
                if self.rows.index(r) == 0 and self.header:
                    html[-1] += '<th %s >%s</th>' %(self.td_flags, cell)
                else:
                    html[-1] += '<td %s >%s</td>' %(self.td_flags, cell)
            html.append(" </tr>")
        html.append("</table>")
        return "\n".join(html) + "\n"

    def asTXT(self):
        """Returns the TEXT representation of the table object."""
        txt = []
        for r in self.rows:
            if self.rows.index(r) == 0 and self.header:
                txt.append('#'+'\t'.join(r))
            else:
                txt.append('\t'.join(r))
        return "\n".join(txt) + "\n"

class orthologyTable(object):
    """Table formatter for orthologs/paralogs."""
    def __init__(self):
        self.rows = []
        self.header = True
        self.table_class = "otable"
        self.header_class = ""
        self.content_class = ""
        self.tr_classes = {}
        self.td_flags = {}
        self.cell_classes = {}
    
    def add_cell(self, row_index, content, classname=None, flag=None):
        """ Add the given object (any text or HTML code) as the content of
        a cell in the row specified by the "row_index" argument. Cells
        are arranged in a stack-like way. """
        # Creates the grid until fitting the given index
        while len(self.rows) <= row_index:
            self.cell_classes[len(self.rows)] = {}
            self.td_flags[len(self.rows)] = {}
            self.rows.append([])
        #add class
        if classname: 
            self.cell_classes[row_index][len(self.rows[row_index])] = classname
        #add flag
        if flag:
            self.td_flags[row_index][len(self.rows[row_index])] = flag
        #add cell content
        self.rows[row_index].append(str(content))

    def add_cells(self, row_index, cells):
        """Add multiple cells."""
        for cell in cells:
            self.add_cell(row_index, cell)

    def asHTML(self):
        """Returns the HTML representation of the table object."""
        html = ['<table class="%s">' % self.table_class]
        tr_class = ""
        for i, row in enumerate(self.rows):
            if i in self.tr_classes:
                tr_class = "class='%s'" % self.tr_classes[i]
            html.append(' <tr %s>\n' % tr_class)
            for j, cell in enumerate(row):
                #get td flag
                if j in self.td_flags[i]:
                    td_flag = self.td_flags[i][j]
                else:
                    td_flag = ""
                #get cell class
                if j in self.cell_classes[i]:
                    cell_class = "class='%s'" %self.cell_classes[i][j]
                else:
                    cell_class = ""                    
                #header
                if   not i and self.header:
                    html[-1] += '<th %s %s>%s</th>' %(cell_class, td_flag, cell)
                else:
                    html[-1] += '<td %s %s>%s</td>' %(cell_class, td_flag, cell)
            html.append(" </tr>")
        html.append("</table>")
        return "\n".join(html) + "\n"

    def asTXT(self):
        """Returns the TEXT representation of the table object."""
        txt = []
        for r in self.rows:
            if self.rows.index(r) == 0 and self.header:
                txt.append('#'+'\t'.join(r))
            else:
                txt.append('\t'.join(r))
        return "\n".join(txt) + "\n"
