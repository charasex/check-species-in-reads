#!/usr/bin/env python3
import sys
import time
import json
import urllib.request
import urllib.parse

def get_ncbi_taxid(species_name, cache={}):
    """
    Interroge l'API NCBI Entrez pour obtenir le TaxID à partir du nom scientifique.
    Utilise un cache local pour éviter les requêtes API répétitives.
    """
    species_clean = species_name.strip()
    if not species_clean or species_clean.lower() in ["sp.", "unknown", "unassigned"]:
        return None

    if species_clean in cache:
        return cache[species_clean]
    
    # Encodage du nom d'espèce pour l'URL
    encoded_name = urllib.parse.quote(species_clean)
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=taxonomy&term={encoded_name}&retmode=json"
    
    try:
        req = urllib.request.Request(url, headers={'User-Agent': 'Kraken2-BOLD-Formatter'})
        with urllib.request.urlopen(req) as response:
            data = json.loads(response.read().decode('utf-8'))
            id_list = data.get('esearchresult', {}).get('idlist', [])
            if id_list:
                taxid = id_list[0]
                cache[species_clean] = taxid
                return taxid
    except Exception as e:
        print(f"[AVERTISSEMENT] Erreur API pour '{species_clean}': {e}", file=sys.stderr)
    
    cache[species_clean] = None
    return None

def format_bold_fasta(input_fasta, output_fasta):
    """
    Formatte un fichier FASTA BOLD pour Kraken 2.
    """
    taxid_cache = {}
    processed_count = 0
    success_count = 0
    missing_species = set()

    with open(input_fasta, 'r') as infile, open(output_fasta, 'w') as outfile:
        current_seq_id = None
        current_species = None
        skip_sequence = False

        for line in infile:
            if line.startswith('>'):
                processed_count += 1
                header_content = line[1:].strip()
                
                # Format standard des entêtes FASTA de BOLD : ProcessID|Species_name|Marker|...
                parts = header_content.split('|')
                
                if len(parts) >= 2:
                    seq_id = parts[0].strip()
                    species_name = parts[1].strip()
                else:
                    # Fallback : prend le 1er mot comme ID et les 2 suivants comme espèce
                    tokens = header_content.split()
                    seq_id = tokens[0]
                    species_name = " ".join(tokens[1:3]) if len(tokens) >= 3 else ""

                # Pause de 0.35s si l'espèce n'est pas encore dans le cache (respect limite 3 req/s NCBI)
                if species_name and species_name not in taxid_cache:
                    time.sleep(0.35)

                taxid = get_ncbi_taxid(species_name, taxid_cache)

                if taxid:
                    new_header = f">{seq_id}|kraken:taxid|{taxid} {species_name}\n"
                    outfile.write(new_header)
                    skip_sequence = False
                    success_count += 1
                else:
                    skip_sequence = True
                    if species_name:
                        missing_species.add(species_name)
                    print(f"[OMIS] TaxID introuvable pour '{species_name}' (ID: {seq_id})", file=sys.stderr)
            else:
                if not skip_sequence:
                    outfile.write(line)

    print("\n================ Bilan du formatage ================")
    print(f"Séquences lues : {processed_count}")
    print(f"Séquences conservées & converties : {success_count}")
    if missing_species:
        print(f"Espèces sans TaxID NCBI ({len(missing_species)}) :")
        for sp in sorted(missing_species):
            print(f"  - {sp}")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python format_bold_kraken.py <input_bold.fasta> <output_kraken.fasta>")
        sys.exit(1)

    format_bold_fasta(sys.argv[1], sys.argv[2])