#/usr/bin/env python
"""Python metaphors db client."""

import MySQLdb, os, re, random, sys
from _mysql_exceptions import OperationalError
from htmlTable import orthologyTable
#import MySQLdb.cursors as cursors
#import word cloud - two paths needed in cgenomics
sys.path.append('/home/lpryszcz/src/python/site-packages')
sys.path.insert(0,'/home/lpryszcz/src/python/site-packages/numpy-1.8.1-py2.6-linux-x86_64.egg')
TAGCLOUD = False
try:
    import wordcloud
    TAGCLOUD = True
except ImportError, e:
    sys.stderr.write("[WARNING] Cannot import wordcloud: %s!\n"%str(e))
#set locale for int printing and so on
import locale
locale.setlocale(locale.LC_ALL, 'en_US.utf8')

class metaphors(object):
    """ Connect to metaPhors database. Return cursor object. """
    def __init__(self, db="metaphors_201601", user="script", passwd="python", \
                 host="cgenomics.crg.es", port=4042):
        self.mysqldb = db
        self.user    = user
        self.passwd  = passwd
        self.host    = host
        self.port    = port
        self.return_html = False
        self.log = sys.stderr
        #function links
        self._fetch = self._execute_and_fetch
        #connect and init db cached objects
        self._connect()
        self._init_db_variables()
        self._init_webserver_vars()
        
    def _connect(self):
        #connect to metaphors
        self.connection = MySQLdb.connect(user=self.user, passwd=self.passwd, \
                                          host=self.host, db=self.mysqldb, port=self.port)
        self.cursor = self.connection.cursor()


    def _execute_and_fetch(self, cmd, fetchOne=False):
        """Execute MySQL cmd and fetch results.
        Reconnect on OperationalError: (2006, 'MySQL server has gone away')
        """
        #execute cmd with exception checking
        try:
            self.cursor.execute(cmd)
        except OperationalError, e:
            #MySQL went away
            if e[0] == 2006:
                self.log.write("Reconnecting to MySQL database...\n")
                self._connect()
                self.cursor.execute(cmd)
            #another error
            else:
                self.log.write('[ERROR] MySQL (%s): %s\n' % (e[0], e[1]))
                return
        #return results
        if fetchOne:
            return self.cursor.fetchone()
        else:
            return self.cursor.fetchall()
        
    def get_spname(self, taxid, short=False):
        """Return species name"""
        spname = self.species[taxid][1]
        if short and len(spname.split())>1:
            spname = "%s. %s"%(spname.split()[0][0], spname.split()[1])
        return spname
            
    def get_random_protid(self):
        """Return random protid"""
        return random.randint(1, self.seqcount + 1)

    def id2int(self, metaid):
        """Return INT metaid"""
        #strip M! and spcode
        if type(metaid) is str:
            #get by extid
            if not metaid.lstrip("M!").split('_')[0].isdigit():
                if metaid.startswith('Phy'): #not self.get_metaid_by_external(metaid):
                    return self.get_metaid_by_external(metaid.split('_')[0])
                return self.get_metaid_by_external(metaid)
            else:
                return int(metaid.lstrip("M!").split('_')[0])
        return metaid

    def id2str(self, metaid, spname=0, spcode=1, uniprot=0, pre="M!"):
        """Return formatted metaid"""
        #strip to int
        if type(metaid) is str:
            metaid = self.id2int(metaid)
        #add sp code
        if spcode:
            taxid = self.get_taxid(metaid)
            code, name = self.species[taxid][:2]
            fullid = "%s%8s_%s"%(pre, metaid, code)
        else:
            fullid = "%s%8s"%(pre, metaid)
        #add gene name and uniprot ids
        if uniprot:
            gene_name = self.get_gene_name(metaid)
            if gene_name:
                fullid += "|%s"%"|".join(gene_name)
            uniprot   = self.get_uniprot_ids(metaid)
            if uniprot:
                fullid += "|%s"%"|".join(uniprot)
        #add sp name
        fullid = fullid.replace(" ","0")
        if spname:
            fullid += " [%s]"%name
        return fullid

    ##GET meta ID
    def get_metaid(self, protid):
        """Return metaid"""
        #clear all blank chars from protid
        if type(protid) is not str:
            return protid
        protid = re.sub(r'\s', '', protid)
        #already metaid given
        if protid.startswith('M!') or protid.isdigit():
            metaid = self.id2int(protid)
        #get by external id
        elif protid.startswith('Phy'): #not metaid and "_" in protid:
            metaid = self.get_metaid_by_external(protid.split('_')[0])
        else:
            metaid = self.get_metaid_by_external(protid)
        return metaid
        
    def get_metaid_by_external(self, protid):
        """Return metaid by externa id"""
        cmd1 = "SELECT protid FROM id_conversion WHERE extid='%s'" % protid
        cmd2 = "SELECT protid FROM ext2meta WHERE extid='%s'" % protid
        for cmd in (cmd1, cmd2):
            if self._fetch(cmd, 1):
                return self._fetch(cmd, 1)[0]
        
    #SEQUENCE RELATED
    def get_species_names(self, term):
        """Return species names with term"""
        return filter(lambda x: term.lower() in x[0].lower(), \
                      sorted(self.name2taxid.iteritems()))
        
    def get_species_info(self, taxid):
        """Return species info"""
        return self._get_genome_info(taxid, True)
                
    def _get_genome_info(self, taxid, more=False):
        """Return stats about genome"""
        cmd = "select count(*) from protid2taxid where taxid=%s"%taxid
        protcount, = self._fetch(cmd, 1)
        out = []
        out.append("%s: <ul><li>%s genes"%(self._add_link(taxid, self.db['taxonomy'][4], \
                                           self.species[taxid][1]), locale.format("%d", protcount, grouping=True)))
        #get more info
        if more:
            #homologs and orthologs
            homologs,  = self._fetch("select count(*) from homologs_%s"%taxid, 1)
            orthologs, = self._fetch("select count(*) from homologs_%s where CS>=0.5"%taxid, 1)
            out.append("%s homologs; %s orthologs (with CS>=0.5)"%(locale.format("%d", homologs, grouping=True), locale.format("%d", orthologs, grouping=True)))
            #external ids
            cmd = """SELECT name, www, version, count(*) FROM ext2taxid
            WHERE taxid=%s GROUP BY name, version ORDER BY dbid, version"""%taxid
            extids = self._fetch(cmd)
            out.append("%s external IDs registered:<ul>"%locale.format("%d", sum(x[-1] for x in extids), grouping=True))
            for dbname, www, version, ids in extids:
                dbname = self._add_link('', www, dbname)
                if version:
                    dbname += " [preoteome version: %s]"%version
                out.append("%s: %s"%(dbname, locale.format("%d", ids, grouping=True)))
            out[-1] += "</ul>"
            
        return "\n<li>".join(out)+"</ul>"

    def get_sequence(self, metaid):
        """Return sequence"""
        cmd="""SELECT seq FROM protid2seq WHERE protid=%s""" % metaid
        return self._fetch(cmd, 1)[0]
        
    def get_taxid(self, metaid):
        """Return taxid associated with protid"""
        cmd = "SELECT taxid FROM protid2taxid WHERE protid=%s"%metaid
        return self._fetch(cmd, 1)[0]

    def get_external_ids(self, metaid):
        """Return external ids"""
        cmd = """SELECT extid, protein_link FROM ext2meta em JOIN db ON em.dbid=db.dbid WHERE protid=%s GROUP BY extid""" % metaid
        return [self._add_link(extid, link) for extid, link in self._fetch(cmd)]

    def get_tree_count(self, metaids):
        """Return number of trees and databases for given metaid."""
        cmd    = "SELECT COUNT(DISTINCT(treeid)), COUNT(DISTINCT(dbid)) FROM protid2treeid WHERE protid=%s"%metaids[0]
        for metaid in metaids[1:]:
            cmd += " AND treeid IN (SELECT treeid FROM protid2treeid WHERE protid=%s)"%metaid
        return self._fetch(cmd, 1)
        
    def get_tree_info(self, metaids):
        """Return formatted info about tree count and db count."""
        #recognise if metaid is list/tuple or single element
        #if isinstance(metaids, basestring) or isinstance(metaids, int) \
        #   or isinstance(metaids, long):
        if not hasattr(metaids, '__iter__'):
            metaids = (metaids,)
        trees, dbs = self.get_tree_count(metaids)
        if not trees:
            return "No trees found!"
        s1 = s2 = ""
        if trees > 1:
           s1 = "s"
        if dbs > 1:
           s2 = "s"
        return "Present in %s phylogenetic tree%s from %s database%s"%(trees, s1, dbs, s2)
        
    def get_uniprot_ids(self, metaid):
        """Return external ids"""
        cmd = """SELECT extid FROM uniprot WHERE protid=%s""" % metaid
        return sorted(set(extid for extid, in self._fetch(cmd)))[:2]
        
    def get_gene_name(self, metaid):
        """Return gene names associated with metaid."""
        cmd = """SELECT extid FROM gene_name WHERE protid=%s""" % metaid
        return sorted(set(extid for extid, in self._fetch(cmd)))[:5]

    def get_gene_names(self, term, limit=100):
        """Return gene names starting with term."""
        cmd = """SELECT name, extid, CONCAT('M!',REPEAT('0',8-LENGTH(protid)),protid,'_',code)
        FROM gene_name_species WHERE extid LIKE '%s%s'
        ORDER BY name, extid LIMIT %s""" % (term, '%', limit)
        return sorted(self._fetch(cmd))

    def get_description(self, metaid):
        """Return gene description"""
        cmd = """SELECT description FROM protid2description WHERE protid=%s"""%metaid
        if self._fetch(cmd, 1):
            return "; ".join(x.split("=")[1].split(';')[0]
                             for x in self._fetch(cmd, 1)[0].split(': ')[1:7:2])
        return 

    def get_GO(self, metaids):
        """Return GO terms associated with metaid(s)."""
        #if type(metaids)==str and ',' in metaids:
        cmd = "SELECT extid, term_type, name FROM go where protid in (%s)"%metaids
        type2go = {}
        for acc, term_type, name in self._fetch(cmd):
            if term_type not in type2go:
                type2go[term_type] = set()
            type2go[term_type].add((name, acc))
        return type2go

    #FORMATTED OUTPUT
    def get_fasta(self, metaid=None):
        """Return FASTA-formatted sequence
        Return RANDOM sequence if no metaid provided."""
        if not metaid:
            metaid = self.get_random_protid()
        metaid = self.id2int(metaid)
        seq    = self.get_sequence(metaid)            
        seqformatted = "\n".join(seq[i:i+60] for i in range(0,len(seq), 60))
        header = self.id2str(metaid, 1, 1, 1)
        return ">%s\n%s" % (header, seqformatted)

    def get_proteome(self, taxid):
        """Return proteme in FASTA"""
        cmd = "select protid from protid2taxid where taxid=%s"%taxid
        return "\n".join(self.get_fasta(metaid) for metaid, in self._fetch(cmd))
        
    def _add_link(self, value, link, value2=""):
        """Add link to property.
        opened in the same window. target='_blank' """
        if not value2:
            value2 = value
        if self.return_html and link:
            if type(value) is str and "|" in value and '%s' in link:
                return "<a href='%s'>%s</a>"%(link%tuple(value.split('|')), value2)
            else:
                return "<a href='%s%s'>%s</a>"%(link, value, value2)
        else:
            return value

    def get_tagcloud(self, metaid, go_terms):
        """Return html link to tag cloud"""
        import matplotlib.pyplot as plt
        htypetxt = ""
        if type(metaid) is str:
            htypetxt = metaid.split("_")[-1]
            metaid = "_".join(metaid.split("_")[:-1])
        metaid = self.id2str(metaid)
        tmpfn = os.path.join(self.TMP, "%s.%s.cloud.png"%(metaid, htypetxt))
        if not os.path.isfile(tmpfn):
            text = []
            for go_data in go_terms.itervalues():
                for name, go in go_data:
                    text.append(name)
            if not text:
                return "No GO terms found!"
            # generate wordcloud
            wc = wordcloud.WordCloud(font_path="DroidSansMono.ttf").generate("\n".join(text))
            wc.to_file(tmpfn)
        return '<img src="%s" title="%s %s functions" width="250px">'%(tmpfn, metaid, htypetxt)
            
    def get_protein_info(self, metaid, short=1):
        """Return protein info"""
        joinby = "\n"
        if self.return_html:
             joinby =  "<br/>"
        metaid  = self.id2int(metaid)
        if not metaid:
            return "Your query was not recognised!"
        fullid      = self.id2str(metaid)
        taxid       = self.get_taxid(metaid)
        gene        = self.get_gene_name(metaid)
        uniprot     = self.get_uniprot_ids(metaid)
        description = self.get_description(metaid)
        #populate
        out = []
        #protein definition - can be obtained with self.id2str(metaid, 1, 1, 1) but without html
        out.append("<table class='protein_info'><tr><td width='750px'><b><a href='/?q=single&metaid=%s'>%s</a></b>"%(fullid, fullid))
        if gene:
            out[-1] += "|<b><i>%s</i></b>" % "|".join(self._add_link("%s|%s"%(g, taxid), \
                                             self.db['gene_name'][4], g) for g in gene)
        if uniprot:
            out[-1] += "|<b>%s</b>" % "|".join(self._add_link(v, self.db['uniprot'][4]) \
                                               for v in uniprot)
        #define organism
        out[-1] += " <b>[%s]</b>"%self._add_link(taxid, self.db['taxonomy'][4], \
                                                 self.species[taxid][1])
        if description:
            out.append("<b>%s</b><br />"%description)
        else:
            out.append("")
        #get tree count
        if short!=4:
            out.append("<b>%s</b>"%self._add_link(fullid, "/?q=tree&metaid=", self.get_tree_info(metaid)))

        #define links
        if short!=3:
            out.append("<b>External identifiers:</b> %s" % ", ".join(self.get_external_ids(metaid)))
                                                    
        out[-1] += "<td>"
        if short in (1, '1'):
            return joinby.join(out)+"</table>"
        #get go terms
        go = []
        go_terms = self.get_GO(metaid)      
        if go_terms:
            if TAGCLOUD:
                out[-1]+=self.get_tagcloud(metaid, go_terms)
            go.append("""<div style='float: left; width: 100%'>
                      <br /><b>Gene Ontology</b><table class='go_table'""")
            html = "<tr><td><b>%s</b><td><ul>%s</ul>\n"
            for term_type, go_data in go_terms.iteritems():
                go.append(html%(term_type.replace("_"," "), 
                                "<li>"+"<li>".join(self._add_link(v, self.db['go'][4], n) \
                                            for n, v in sorted(go_data))))
            go.append("</table></div>")
        out[-1]+="<div id='second_wordcloud'></div></table>"
        #return without sequence
        if   short>1: 
            return joinby.join(out) 
        #define sequence
        out.append("<div style='float: left'><b>Sequence in FASTA:</b>")
        fasta   = self.get_fasta(metaid)
        if self.return_html:
             fasta  = "<pre>%s</pre></div>" % fasta
        return joinby.join(out) + "\n" + fasta + "\n" + "\n".join(go)
        
    def get_trees(self, ids):
        """Return trees for protein IDs."""
        out = []
        if type(ids) is str:
            ids = ids.split()
        #get internal ids
        metaids = [self.id2int(x) for x in ids if x and self.id2int(x)]
        if not metaids:
            return "Your query was not recognised!"
        cmd = """SELECT name, www, protein_link, tree_link, extid, treeid
                 FROM protid2treeid WHERE protid=%s"""%metaids[0]
        for i, metaid in enumerate(metaids):
            out.append(self.get_protein_info(metaid, 3))
            #limit by common with other metaids
            if i:
                cmd += " AND treeid IN (SELECT treeid FROM protid2treeid WHERE protid=%s)"%metaid
        #remove repeats
        cmd += " GROUP BY treeid"
        results = self._fetch(cmd)
        ##format output
        if len(metaids)>1:
            text = "<b><table><tr><td>Selected proteins are together %s</table></b>"
            out.append(text%self.get_tree_info(metaids).lower())
        if not results:
            return "".join(out)
        #process
        db2trees = {}
        for dbname, www, protein_link, tree_link, extid, treeid in results:
            dblink = self._add_link('', www, dbname)
            #add db
            if dblink not in db2trees:
                db2trees[dblink] = {}
            #add extid
            extid = self._add_link(extid, protein_link)
            if extid not in db2trees[dblink]:
                db2trees[dblink][extid] = []
            #add trees
            db2trees[dblink][extid].append(self._add_link(treeid, tree_link))
        #get html
        table = orthologyTable()
        #add header
        colnames = ("Database", "No. of trees", "External ID", "Trees")
        widths = (80, 80, 120, 0) 
        for cell, width in zip(colnames, widths):
            cell_class = td_flag = ""
            if width:
              td_flag = 'width="%s"'%width
            table.add_cell(0, cell, cell_class, td_flag)
        row_index = 1
        for i, db in enumerate(sorted(db2trees), 1):
            td_flag = 'rowspan="%s"' % len(db2trees[db])            
            #even or odd species row
            if i%2:
                table.tr_classes[i] = "even_row"
            else:
                table.tr_classes[i] = "odd_row"
            #add db & tree count cell
            cell_class = None
            table.add_cell(row_index, db, cell_class, td_flag)
            table.add_cell(row_index, sum(len(x) for x in db2trees[db].itervalues()), cell_class, td_flag)
            #add extid & trees
            for extid, trees in sorted(db2trees[db].iteritems()):
                table.add_cell(row_index, extid)
                table.add_cell(row_index, "<div id='scrollcell'>"+", ".join(sorted(trees))+"</div>")
                row_index += 1

        if self.return_html:
            return "<br />".join(out) + table.asHTML()
        else:
            return "\n#".join(out) + table.asTXT()

    def get_orthologs_table(self, protid):
        """Return HTML or TXT formatted orthologs table"""
        htypetxt = "ortholog"
        #get metaid
        metaid = self.get_metaid(protid)
        ##catch wrong protein ids here
        if not metaid:
            return "Your query produced no hits!"
        #add short protein info
        header = self.get_protein_info(metaid, 2)
        #get orthologs
        orthologs = self.get_orthologs(metaid)
        if not orthologs:
            return header+"<b>No %ss have been found.</b>"%htypetxt
        header += "<b>%s %ss in %s species</b>\n"%(sum(len(o) for o in orthologs.itervalues()), \
                                                   htypetxt, len(orthologs))
        out = []
        otable = orthologyTable()
        #add header
        colnames = ("Target species", "Co-%ss"%htypetxt, "%ss"%htypetxt.capitalize(), \
                    "CS", "EL", "Trees")+ self.TREES_REPOSITORIES_NAMES
        widths = (180, 120, 120, 0, 0, 0) + (0,)*len(self.TREES_REPOSITORIES_NAMES)
        for cell, width in zip(colnames, widths):
            cell_class = td_flag = ""
            if width:
                td_flag = 'width="%s"'%width
            otable.add_cell(0, cell, cell_class, td_flag)
        #populate table
        row_index = 1
        taxids = sorted([t for t in orthologs if t in self.species], key=lambda x: self.species[x][1])
        for i, taxid in enumerate(taxids, 1):
            td_flag = 'rowspan="%s"' % len(orthologs[taxid])
            spname = self.species[taxid][1]
            cell_class = None
            otable.add_cell(row_index, spname, cell_class, td_flag)
            #even or odd species row
            if i%2:
                otable.tr_classes[row_index] = "even_row"
            else:
                otable.tr_classes[row_index] = "odd_row"
            #add orthologs
            for j, (metaid2, extid2, cs, el, trees, css, coorths) in enumerate(orthologs[taxid], 1):
                #use metaid if not external id
                if not extid2:
                    extid2 = self.id2str(metaid2, 0, 0)
                #co-orthology information
                coorthologs = self.id2str(metaid, 0, 0) + " " + \
                              " ".join(self.id2str(coid, 0, 0) for coid, cocs in coorths)
                #collapse co-orthologs if too many
                if self.return_html and coorths:
                    coorthologs = """
   <a onClick="javascript:togglecomments(\'i%s_j%s\')"> %s co-%ss</a>
    <div id="i%s_j%s" class="commenthidden">%s</div>
   </a>""" % (i, j, len(coorths)+1, htypetxt, i, j, coorthologs)
                otable.add_cell(row_index, coorthologs)
                #ortholog - popup on scroll over
                if self.return_html:
                    #protid2 = self._add_link(extid2, self.db['metaphors'][4])
                    protid2 = self.SEQLINKOR % (metaid2, metaid2, extid2, metaid2)
                    trees = self._add_link("%s|%s"%(self.id2str(metaid), self.id2str(metaid2)), \
                                           self.db['metaphors'][3], trees)
                else:
                    protid2 = self.id2str(metaid2)
                #protid2, CS,  EL & trees
                otable.add_cells(row_index, (protid2, "%.3f"%cs, el, trees))
                #tree repository signals
                for signal in css:
                    cell = cell_class = ""
                    if signal:
                        extcs, exttrees = signal
                        if exttrees>1:
                            cell = "%.2f/%s" % (extcs, exttrees)
                        if extcs >= self.csth:
                            cell_class = "agreement"
                        else:
                            cell_class = "disagreement"
                    otable.add_cell(row_index, cell, cell_class)
                row_index += 1
        #return html-formatted
        if self.return_html:
            return header + otable.asHTML() + self.LEGEND
        #OrthoXML to ADD
        #TXT
        else:
            return otable.asTXT()

    ##GET ORTHOLOGS
    def get_paralogs(self, protid, csth=0.5, elth=1, extdbs=[], homtype=0, compare="<"):
        """Return dictionary of paralogs.

        Each dictionary contain taxid as key and then list of orthologs/paralogs
        is given. For each ortholog/paralog, metaid, CS,
        database signals (dbid:orhtology_trees,all_trees) and list of
        co-orthologs/co-paralogs with respective consistecy scores are given.
        """
        return self.get_orthologs(protid, 1-csth, elth, extdbs, homtype, compare)
        
    def get_orthologs(self, protid, csth=0.5, elth=1, extdbs=[], homtype=1, compare=">="):
        """Return dictionary of orthologs.
        
        Each dictionary contain taxid as key and then list of orthologs/paralogs
        is given. For each ortholog/paralog, metaid, CS,
        database signals (dbid:orhtology_trees,all_trees) and list of
        co-orthologs/co-paralogs with respective consistecy scores are given.
        """
        homologs = {}
        #get taxid and internal query
        protid1 = self.get_metaid(protid)
        if not protid1: 
            return homologs
        taxid = self.get_taxid(protid1)
        cmd = """SELECT taxid, h.protid2, i.extid, h.CS, h.sources, h2.protid1, h2.CS
        FROM homologs_%s AS h JOIN homologs_%s AS h2 ON h.protid2=h2.protid2
        JOIN protid2taxid AS t ON h.protid2=t.protid
        LEFT JOIN uniprot AS i ON h.protid2=i.protid
        WHERE h.protid1 = %s AND h.CS %s %s AND h2.CS %s %s
        """%(taxid, taxid, protid1, compare, csth, compare, csth)
        pprotid2 = [] 
        for taxid2, protid2, extid2, cs, sources, coprotid, cocs in self._fetch(cmd): 
            #skip results not passing EL criteria etc
            include = self._include(sources, csth, elth, extdbs, homtype)
            if not include: 
                continue
            trees, el, css = include
            #add homologs info only once for each protid2
            if not pprotid2 or protid2 != pprotid2[-1]:
                #add taxa if not yet added
                if taxid2 not in homologs:
                    homologs[taxid2] = []
                #store in homologs
                homologs[taxid2].append([protid2, extid2, cs, el, trees, css, []])
                #update pprotid
                pprotid2.append(protid2)
            #add co-homologs different than query! # show_structure == '1' 
            if coprotid != protid1:
                homologs[taxid2][-1][-1].append((coprotid, cocs))
        return homologs

    def get_orthologs_and_paralogs(self, protid, csth=0.5, elth=1, extdbs=[]):
        """Return dictionary of orthologs and paralogs.
        Each dictionary contain taxid as key and then list of orthologs/paralogs
        is given. For each ortholog/paralog, metaid, CS,
        database signals (dbid:orhtology_trees,all_trees) and list of
        co-orthologs/co-paralogs with respective consistecy scores are given.
        """
        orthologs, paralogs = {}, {}
        #get taxid and internal query
        protid1 = self.get_metaid(protid)
        if not protid1: 
            return orthologs, paralogs
        taxid = self.get_taxid(protid1)
        cmd = """SELECT taxid, h.protid2, i.extid, h.CS, h.sources, h2.protid1, h2.CS
        FROM homologs_%s AS h JOIN homologs_%s AS h2 ON h.protid2=h2.protid2
        JOIN protid2taxid AS t ON h.protid2=t.protid
        LEFT JOIN uniprot AS i ON h.protid2=i.protid
        WHERE h.protid1 = %s""" % (taxid, taxid, protid1)
        pOprotid2, pPprotid2 = [], []
        for taxid2, protid2, extid2, cs, sources, coprotid, cocs in self._fetch(cmd): 
            #use orthologs/paralogs and pprotid from orthologs
            if   cs >= csth and cocs >= csth:
                homologs = orthologs
                pprotid2 = pOprotid2
                homtype  = 1
            #or paralogs depending on cs
            elif cs < csth and cocs < csth:
                homologs = paralogs
                pprotid2 = pPprotid2
                homtype  = 0
            else:
                continue
            #skip results not passing EL criteria etc
            include = self._include(sources, csth, elth, extdbs, homtype)
            if not include: 
                continue
            trees, el, css = include
            #add homologs info only once for each protid2
            if not pprotid2 or protid2 != pprotid2[-1]:
                #add taxa if not yet added
                if taxid2 not in homologs:
                    homologs[taxid2] = []
                #store in homologs
                homologs[taxid2].append([protid2, extid2, cs, el, trees, css, []])
                #update pprotid
                pprotid2.append(protid2)
            #add co-homologs different than query! # show_structure == '1' 
            if coprotid != protid1:
                homologs[taxid2][-1][-1].append((coprotid, cocs))
        return orthologs, paralogs 

    def get_orthologs_between_two_species_table(self, taxid1, taxid2, csth=0.5, elth=1):
        """Return formatted orthology table between two species """
        # 1.3M entries for YEAST-vs-CANAL due to id_conversion                                                                           
        return "Work in progress... Coming soon!"          
        #convert species names into taxid
        if not str(taxid1).isdigit():
            if not taxid1 in self.name2taxid:
                return "Species %s have not been found!"%taxid1
            taxid1 = self.name2taxid[taxid1]
        if not str(taxid2).isdigit():
            if not taxid2 in self.name2taxid:
                return "Species %s have not been found!"%taxid2
            taxid2 = self.name2taxid[taxid2]
        htypetxt = "ortholog"
        #add short species info ADD protein counts etc
        header = self._get_genome_info(taxid1)
        header += self._get_genome_info(taxid2)        
        #get orthologs
        orthologs = self.get_orthologs_between_two_species(taxid1, taxid2)
        ##THIS IS MESSY BUT WORKING!
        ocount  = sum(len(o) for o in orthologs.itervalues())
        p2count = {}
        for o in filter(lambda x: len(x)==1, orthologs.itervalues()):
            p = o[0][0]
            if p not in p2count:
                p2count[p]  = 1
            else:
                p2count[p] += 1
        for o in filter(lambda x: len(x)>1, orthologs.itervalues()):
            for oi in o:
                p = oi[0]
                if p in p2count:
                    p2count[p] += 1            
        ocount_one2one = sum(c for p, c in p2count.iteritems() if c==1)
        header += "<b> %s orthologs of which %s are one-to-one.</b>"%(locale.format("%d", ocount, grouping=True), locale.format("%d", ocount_one2one, grouping=True))
        #get orthologs
        otable = orthologyTable()
        #add header
        colnames = ("<i>%s</i>"%self.get_spname(taxid1, 1), \
                    "<i>%s</i>"%self.get_spname(taxid2, 1), \
                    "CS", "EL", "Trees")+ self.TREES_REPOSITORIES_NAMES
        widths = (120, 120, 0, 0, 0) + (0,)*len(self.TREES_REPOSITORIES_NAMES)
        for cell, width in zip(colnames, widths):
            cell_class = td_flag = ""
            if width:
                td_flag = 'width="%s"'%width
            otable.add_cell(0, cell, cell_class, td_flag)
        #populate table
        row_index = 1
        for i, (prot1data, odata) in enumerate(orthologs.iteritems()):
            #(metaid1, extid1)
            metaid1 = int(prot1data.split('_')[0])
            extid1  = "_".join(prot1data.split('_')[1:])
            td_flag = 'rowspan="%s"' % len(odata)
            cell_class = None
            #use metaid if not external id
            if not extid1 or extid1=="None":
                extid1 = self.id2str(metaid1, 0, 0)
            #ortholog - popup on scroll over
            if self.return_html:
                protid1 = self.SEQLINKOR % (metaid1, metaid1, extid1, metaid1)
            else:
                protid1 = self.id2str(metaid1)
            otable.add_cell(row_index, protid1, cell_class, td_flag)
            #even or odd species row
            if i%2:
                otable.tr_classes[row_index] = "even_row"
            else:
                otable.tr_classes[row_index] = "odd_row"
            for metaid2, extid2, cs, el, trees, css in odata:
                #use metaid if not external id
                if not extid2:
                    extid2 = self.id2str(metaid2, 0, 0)
                #ortholog - popup on scroll over
                if self.return_html:
                    protid2 = self.SEQLINKOR % (metaid2, metaid2, extid2, metaid2)
                    trees = self._add_link("%s|%s"%(self.id2str(metaid1), self.id2str(metaid2)), \
                                           self.db['metaphors'][3], trees)
                else:
                    protid2 = self.id2str(metaid2)
                #protid2, CS,  EL & trees
                otable.add_cells(row_index, (protid2, "%.3f"%cs, el, trees))
                #tree repository signals
                for signal in css:
                    cell = cell_class = ""
                    if signal:
                        extcs, exttrees = signal
                        if exttrees>1:
                            cell = "%.2f/%s" % (extcs, exttrees)
                        if extcs >= self.csth:
                            cell_class = "agreement"
                        else:
                            cell_class = "disagreement"
                    otable.add_cell(row_index, cell, cell_class)
                row_index += 1
        #return html-formatted
        if self.return_html:
            return header + otable.asHTML() + self.LEGEND
        #OrthoXML to ADD
        #TXT
        else:
            return otable.asTXT()               
            
    ##ORTHOLOGS/PARALOGS FOR TWO SPECIES 
    def get_paralogs_between_two_species(self, taxid1, taxid2, csth=0.5, elth=1, \
                                         extdbs=[]):
        """Return paralogs between two taxa. Alias function."""
        return self.get_orthologs_between_two_species(taxid1, taxid2, csth, elth, \
                                                      extdb, homtype=0)
        
    def get_orthologs_between_two_species(self, taxid1, taxid2, csth=0.5, elth=1,\
                                          extdbs=[], homtype=1):
        """Return all orthologs between two taxa.
        Return paralogs if homtype == 0
        """
        orthologs, p2coorthologs = {}, {}
        compare = '>='
        if homtype == 0:
            compare = '<' #searching for paralogs
        cmd="""SELECT protid1, p1.extid, protid2, p2.extid, CS, sources
        FROM `homologs_%s` AS h JOIN protid2taxid AS t ON h.protid2=t.protid
        LEFT JOIN uniprot AS p1 ON protid1=p1.protid
        LEFT JOIN uniprot AS p2 ON protid2=p2.protid
        WHERE t.taxid=%s AND CS %s %s""" % (taxid1, taxid2, compare, csth)
        for protid1, extid1, protid2, extid2, cs, sources in self._fetch(cmd):
            protid1 = "%s_%s"%(protid1, extid1)
            #skip results not passing EL criteria etc
            include = self._include(sources, csth, elth, extdbs, homtype)
            if not include: 
                continue
            trees, el, css = include
            data_list = [protid2, extid2, cs, el, trees, css]#, []]
            #store info
            if not protid1 in orthologs:
                orthologs[protid1] = []
            orthologs[protid1].append(data_list)
        return orthologs
        
    def _include(self, sources, csth, elth, extdbs, homtype, exclude_db_list=[]):
        """Return true if signals passed by external_data_list are
        following criteria in exclude_db_list and evidence_level
        """
        _include = False
        el = tottrees = 0
        external_data_list=[[] for i in range(len(self.TREES_REPOSITORIES))]
        #8:0.818,22;10:1,1
        for source in sources.split(';'):
            dbid, pred = source.split(':')
            dbid  = int(dbid)
            cs, t = pred.split(',')
            if not cs:
                cs = 0
            cs, t = float(cs), int(t)
            #get position in external data list
            dbcode = self.id2db[dbid]
            #ensembl_* -> ensembl
            dbcode = dbcode.split('_')[0] 
            idx = self.TREES_REPOSITORIES.index(dbcode)
            #add info
            external_data_list[idx].append((cs, t))
            tottrees += t
            
        #collapse multiple sources for one db (ensembl)
        for i, sData in enumerate(external_data_list):
            if not sData:
                continue
            if len(sData)>1:
                otrees = trees = 0
                for cs, t in sData:
                    otrees += round(cs/t)
                    trees  += t
                cs = otrees * 1.0 / trees
            else:
                cs, trees = external_data_list[i][0]
            # paralogy
            if homtype == 0:
                cs = 1-cs
            external_data_list[i] = (cs, trees)

            
        #CHECK IF GIVEN DATABASE CONFIRM PREDICTION
        i = 0
        for dbcode in self.TREES_REPOSITORIES:
            if dbcode in extdbs:
                signal = external_data_list[i]
                if not signal:
                    return
                cs, trees = signal
                if cs < csth:
                    # reject prediction if it's not confirmed by particular db, or phylomes no < 1
                    return
            i += 1
            
        # CHECK EVIDENCE LEVEL AND exclude_db_list criterium
        css = []
        for dbcode, signal in zip(self.TREES_REPOSITORIES, external_data_list):
            if not signal:
                css.append(-1)
                continue
            cs, trees = signal
            css.append(cs)
            if dbcode not in exclude_db_list:
                if cs > csth:
                    _include = True
                    el += 1
        
        #DON'T INCLUDE IF EL_TH > that no of dbs giving positive signals
        if el < elth:
            return

        if _include:
            return tottrees, el, external_data_list

    ###VARIABLES SETTINGS
    def _init_db_variables(self):
        """Load db info, species etc"""
        self._load_db_info()
        self._load_species_info()
        #get sequence count
        self.seqcount, = self._fetch("select count(*) from protid2seq", 1)
        #get protein db length for E-value calculation in rapsi
        ##factor of 45.4... need to be adjusted if table schema changes
        cmd1 = "SELECT DATA_LENGTH-TABLE_ROWS*45.499855814 FROM information_schema.TABLES WHERE TABLE_NAME='protid2seq'"
        self.dblength, = self._fetch(cmd1, 1)
        self.dblength  = int(round(self.dblength))
        
    def _load_db_info(self):
        """Initialise dictionary of all db registered in metaphors.""" 
        self.db, self.dbid2info = {}, {}
        self.db2id = {} #deprecated, but needed for back-compatibility
        self.id2db = {}
        cmd = "SELECT dbid, code, name, www, tree_link, protein_link FROM db"
        for dbid, code, name, www, tree_link, protein_link in self._fetch(cmd):
            self.dbid2info[dbid] = (code, name, www, tree_link, protein_link)
            self.db[code] = (dbid, name, www, tree_link, protein_link)
            self.db2id[code] = dbid
            self.id2db[dbid] = code
            
    def _load_species_info(self):
        """Initialise dictionary of all species available in db.
        taxid2info = { taxid1: ( scientific name,common name,lineage ), ...""" 
        self.species = {}
        self.name2taxid = {}
        cmd = "SELECT taxid, code, name, common, lineage FROM species"
        for taxid, code, name, common, lineage in self._fetch(cmd):
            self.species[taxid] = (code, name, taxid, common, lineage)
            self.name2taxid[name] = taxid
    
    def _init_webserver_vars(self):
        """Set variables, previously stored in _vars.py"""
        #default session
        #NEED TO READ FROM cookies HERE
        self.csth  = 0.5
        self.elth  = 1
        self.paralogs = False
        #global vars        
        self.BIN, self.TMP = "bin", "tmp"
        #similaritySearchNew
        self.BLATBIN = os.path.join(self.BIN, 'blat')
        self.KMER, self.STEP, self.SEQLIMIT, self.HTMLOUT = 5, 1, 100, 1,
        self.SAMPLING = (50, 200)
        self.HASHTABLE = "hash2protids"
        #self.SEQCMD    = "select concat(p.protid,'_',code), seq from protid2taxid p join protid2seq ps join species s on p.protid=ps.protid and p.taxid=s.taxid where p.protid in (%s)"
        self.SEQCMD    = "SELECT CONCAT(REPEAT('0',8-LENGTH(ps.protid)),ps.protid,'_',code) AS protid, seq FROM protid2taxid p JOIN protid2seq ps ON p.protid=ps.protid JOIN species s ON p.taxid=s.taxid WHERE p.protid in (%s)"
        self.SEQLINK   = """<a id="%s" onmouseover="ShowPopUp(this)" onmouseout="ShowPopUp(this)" href="/?q=single&metaid=M!%s">M!%s</a><div id="%s_popUp" class="popUp"></div>"""
        self.SEQLINKOR = """<a id="%s" onmouseover="ShowPopUp(this)" onmouseout="ShowPopUp(this)" href="/?q=sequence&metaid=M!%s">%s</a><div id="%s_popUp" class="popUp"></div>"""
        #orthology table - NEED REDESIGN!
        self.TREES_REPOSITORIES = ('phylome', 'ensembl', 'eggnog',  'orthomcl', 'treefam', 'hogenom')
        self.TREES_REPOSITORIES_NAMES = ('Phylome', 'Ensembl', 'EggNOG', 'OrthoMCL', 'TreeFAM', 'Hogenom')
        self.TREES_REPOSITORIES_FULLNAMES = ('PhylomeDB', 'Ensembl', 'EggNOG', 'OrthoMCL', 'TreeFAM', 'Hogenome')
        self.EXTERNAL_TREES_REPOSITORIES = self.TREES_REPOSITORIES[1:]
        self.EXTERNAL_REPOSITORIES_NAMES = self.TREES_REPOSITORIES_NAMES[1:]
        self.LEGEND = """\n<div id='legend'>
Above table(s) list predicted orthologs and/or paralogs. Predictions were retrieved
from PhylomeDB and other homology repositories. Consistency score, number of trees,
and signals from each repository are given for each prediction - green for consistent
prediction, red for inconsistencies; Ens - Ensembl, Egg - EggNOG, Ort - OrthoMCL,
Hg - Hogenom, TF - TreeFAM.
In addition, co-orthologs from query species are given. On this basis,
one-to-one, one-to-many, and many-to-many relationships can be distinguished.
<br />For more information click <a target=_blank href='http://orthology.phylomedb.org/?q=help'>here</a>.
</div>
"""
        self.VERBOSE = 1
    
    ###OPTIONS LAYER
    '''def get_options(self):
        """ """
        html = []
        joinby="\n"
        if self.return_html:
            joinby="\n<br />"
        for k, v in zip(("CS cut-off", "EL cut-off",), \
                        (self.csth, self.elth)):
            html.append("%s: %s"%(k, v))
        return joinby.join(html)
    '''
    def options(self, cookie):
        """Handle session options."""
        show_structure_text=o_text=p_text=''
        if cookie['show_structure'].value=='true': show_structure_text=' checked'
        if cookie['homology_type'].value=='1':     o_text=' selected'
        else: p_text=' selected'

        html_text=txt_text=xml_text=''
        if   cookie['format'].value=='html': html_text=' selected'
        elif cookie['format'].value=='txt':  txt_text=' selected'
        elif cookie['format'].value=='xml':  xml_text=' selected'

        out="""
      <div id="left_pannel style="display: inline">
        NOT FINISHED!<br />
                          Get:
                          <select size="1" name="homology_type">
                            <option value="1"%s>orthologs</option>
                            <option value="0"%s>paralogs</option>
                          </select> <br />
                          CS cut-off: <input name="CS_th" type="text" value="%s" size="2"> <br />
                          Evidence level: <input name="EL_th" type="text" value="%s" size="2"><br />
                          <input type="checkbox" name="show_structure" value="1"%s> Show homology structure
      </div>
      <div id="right_pannel" style="display: block">
      Predictions confirmed by:<br />	""" % ( o_text,p_text,cookie['CS_th'].value, cookie['EL_th'].value,show_structure_text )

        for db, db_name in zip(self.TREES_REPOSITORIES, self.TREES_REPOSITORIES_FULLNAMES):
            db_text=''
            if db in cookie['external_dbs'].value: db_text=' checked'
            out+="""
                        <input type="checkbox" name="external_dbs" value="%s"%s> %s<br />	""" % ( db,db_text,db_name )


        out+="""

                        <b>Retrieve as:</b>
                          <select size="1" name="format">
                            <option value="html"%s>html</option>
                            <option value="txt"%s>text</option>
                            <!-- <option value="xml"%s>text</option>(not implemented yet) -->
                          </select>
        </div>
        <button name="buttonSubmit" onClick="options('click');">Save</button>
        <button name="buttonReset" onClick="reset_session('%s');">Reset session</button>

                          """ % (html_text, txt_text, xml_text, cookie['sid'].value)

        return out
    