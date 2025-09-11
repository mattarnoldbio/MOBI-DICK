
import pandas as pd
from bs4 import BeautifulSoup
import argparse
import os

def find_taxon(hit, search_rank):
    rank = ""
    taxon = ""
    try:
        while rank != search_rank:
            hit = hit.findParent()
            rank = hit.find("rank").find("val").text
        taxon = hit.attrs["name"]
    except:
        taxon = "none"
    return taxon


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Process Krona Plot to CSV')
    parser.add_argument("-k","--krona_plot", help='Path to Krona plot', default=None)
    parser.add_argument("-c","--contig_file", help='Path to contig file (fasta format)', default=None)
    parser.add_argument("-o","--output_dir", help='Directory to save output to', default=None)
    parser.add_argument("-s","--score_filter", help='Filter hits by mean score', default=-10.0, type=float)
    parser.add_argument("-w","--which_db", help='Which database was used for the Krona plot', default="diamond")
    parser.add_argument("-m","--krona_metadata_dir", help='Path to directory containing krona metadata (if too many hits were generated, this will have been created by kronatools)', default=None)
    args = parser.parse_args()

    krona_plot_file = args.krona_plot
    contig_file = args.contig_file
    output_dir = args.output_dir
    which_db = args.which_db
    krona_metadata_dir = args.krona_metadata_dir
    accession = contig_file.split("/")[-1].split(".")[0]
    out_filename = accession
    if krona_metadata_dir is None:
        krona_metadata_dir = output_dir + "/" + accession + "." + which_db + ".krona.html.files/" 
    
    score_filter = float(args.score_filter)
    if output_dir == "":
        output_dir = contig_file.strip(contig_file.split("/")[-1])

    with open(contig_file) as fp:
        contigs = fp.readlines()

    ids = {}

    for line in contigs:
        if line[0] == ">":
            line_=line.strip(">").strip("\n").split(" ")
            line__ = [ float(x.split("=")[1]) for x in line_[1:]][1:]
            ids[line_[0]] = line__


    with open(krona_plot_file) as fp:
        soup = BeautifulSoup(fp, "html.parser")

    virus_hits = []

    try:
        virus_hits = soup.find_all(attrs={"name": "Viruses"})[0].find_all("members")
        print("Found {} virus hits for {}".format(len(virus_hits), accession))


    except:
        os.system("touch {}/{}_{}_no_virus_hits.txt".format(output_dir, which_db, accession))
        print("No virus hits found for {}".format(accession))
        if which_db != "blastn":
            exit(1)

    if len(virus_hits) > 0:
        hits_df = pd.DataFrame(columns=["sample","contig","species", "genus", "family", "mean_score"])

        for i, virus in enumerate(virus_hits):
            species = virus.findParents()[0].attrs["name"]
            #print(virus.findParents().find_all(attrs={"rank" : "genus"}))
            genus = find_taxon(virus, "genus")
            family = find_taxon(virus, "family")
            mean_score = virus.findParents()[0].find("score").text
            for hit in virus.findParents()[0].find_all("members")[0].find_all("val"):
                contig = hit.text
                if "|" in contig: 
                    accession = contig.split("|")[0]   
                    contig_id = contig.split("|")[1]
                else:
                    contig_id = contig
                try: 
                    #ids[contig_id]
                    hits_df.loc[len(hits_df)] = [accession, contig_id, species, genus, family, float(mean_score)]
                except:
                    contig_metadata = krona_metadata_dir + "/" + contig
                    with open(contig_metadata) as fp:
                        lines = fp.readlines()
                    lines[0] = lines[0].strip("data('")
                    lines = [x.replace(r"\n","") for x in lines]
                    lines = [x.replace("\\","") for x in lines]
                    lines = [x.replace("\n","") for x in lines][:-1]
                    for line in lines:
                        #lsprint(contig)
                        contig = line
                        if "|" in contig:    
                            accession = contig.split("|")[0]  
                            contig_id = contig.split("|")[1]
                        else:
                            contig_id = contig
                        hits_df.loc[len(hits_df)] = [accession, contig_id, species, genus, family, float(mean_score)]


        hits_df.loc[hits_df["mean_score"] < score_filter].to_csv("{}/{}_{}_virus_hits.csv".format(output_dir, which_db, out_filename), index=False)



    if which_db == "blastn":
        all_hits=soup.find_all("members")
        non_viral_hits = []
        for hit in all_hits:
            if find_taxon(hit, "superkingdom") != "Viruses":  
                non_viral_hits.append(hit)

        non_viral_hits_df = pd.DataFrame(columns=["sample","contig","species", "genus", "family", "mean_score"])


        for i, virus in enumerate(non_viral_hits):
            species = virus.findParents()[0].attrs["name"]
            #print(virus.findParents().find_all(attrs={"rank" : "genus"}))
            genus = find_taxon(virus, "genus")
            family = find_taxon(virus, "family")
            if species != "Root":
                mean_score = virus.findParents()[0].find("score").text
            else:
                mean_score = "0"            
            for hit in virus.findParents()[0].find_all("members")[0].find_all("val"):
                contig = hit.text
                if "|" in contig:    
                    accession = contig.split("|")[0]  
                    contig_id = contig.split("|")[1]
                else:
                    contig_id = contig
                try: 
                    non_viral_hits_df.loc[len(non_viral_hits_df)] = [accession, contig_id, species, genus, family, float(mean_score)]
                except:
                    contig_metadata = krona_metadata_dir + "/" + contig
                    with open(contig_metadata) as fp:
                        lines = fp.readlines()
                    lines[0] = lines[0].strip("data('")
                    lines = [x.replace(r"\n","") for x in lines]
                    lines = [x.replace("\\","") for x in lines]
                    lines = [x.replace("\n","") for x in lines][:-1]
                    for line in lines:
                        #print(contig)
                        contig = line
                        if "|" in contig:    
                            accession = contig.split("|")[0]  
                            contig_id = contig.split("|")[1]
                        else:
                            contig_id = contig
                        non_viral_hits_df.loc[len(non_viral_hits_df)] = [accession, contig_id, species, genus, family, float(mean_score)]

        non_viral_hits_df.loc[non_viral_hits_df["mean_score"] < score_filter].to_csv("{}/{}_{}_non_virus_hits.csv".format(output_dir, which_db, out_filename), index=False)
