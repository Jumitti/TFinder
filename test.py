import requests
import time

def get_gene_info(gene_id):
    """
    Récupère les informations d'un gène via l'API NCBI Entrez en utilisant l'ID du gène.
    """
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
    params = {
        "db": "gene",  # On cherche dans la base "gene"
        "id": gene_id,  # ID du gène obtenu via esearch
        "retmode": "json"  # Retourne le résultat en format JSON
    }

    response = requests.get(url, params=params)

    # Vérifier si la requête est réussie
    if response.status_code == 200:
        print(response.text)
        return response.json()
    else:
        print(f"Erreur : Impossible de récupérer les informations (status code: {response.status_code})")
        return None


def extract_genomic_info(gene_info):
    print(gene_info)
    time.sleep(1)
    """
    Extrait les NC_, chrstart et chrstop de la réponse JSON.
    Ne garde que les accession versions qui partagent la même base avant le point.
    """
    accession_dict = {}

    # Extraction depuis la section locationhist
    location_hist = gene_info.get('locationhist', [])
    if len(location_hist) == 0:
        location_hist = gene_info.get('genomicinfo', [])
    for loc in location_hist:
        nc_accver = loc.get('chraccver')  # Extrait le NC_XXXXX.YY
        chrstart = loc.get('chrstart')
        chrstop = loc.get('chrstop')
        print(nc_accver, chrstart, chrstop)

        if nc_accver:
            base_accession = nc_accver
            if base_accession not in accession_dict:
                accession_dict[base_accession] = (chrstart, chrstop)
            else:
                # Conserver la première occurrence des coordonnées
                existing_start, existing_stop = accession_dict[base_accession]
                accession_dict[base_accession] = (min(existing_start, chrstart), max(existing_stop, chrstop))

    print(accession_dict)  # Afficher le dictionnaire des NC_ sans les suffixes
    return accession_dict


def fetch_nc_info(nc_accver):
    time.sleep(1)
    """
    Récupère les informations détaillées pour un NC_ via l'API NCBI.
    """
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
    params = {
        "db": "nuccore",
        "id": nc_accver,
        "retmode": "json"
    }

    response = requests.get(url, params=params)

    if response.status_code == 200:
        print(response.text)
        return response.json()
    else:
        print(f"Erreur : Impossible de récupérer les informations pour {nc_accver} (status code: {response.status_code})")
        return None


def main():
    gene_id = input("Entrez l'ID du gène (par ex. 4843 pour NOS2) : ").strip()

    # Étape 1 : Récupérer les informations du gène via son ID
    gene_info = get_gene_info(gene_id)

    if gene_info and 'result' in gene_info and gene_id in gene_info['result']:
        gene_details = gene_info['result'][gene_id]
        print(gene_details)
        print(f"Informations pour le gène (ID: {gene_id}):")

        # Étape 2 : Extraire les NC_, chrstart et chrstop
        nc_dict = extract_genomic_info(gene_details)
        print("\nNC_ Accession, chrstart, chrstop:")

        # Affichage des NC_ et leurs informations de position
        for base_accver, (chrstart, chrstop) in nc_dict.items():
            print(f"{base_accver}: chrstart={chrstart}, chrstop={chrstop}")

        # Filtrer pour garder uniquement les entrées qui commencent par "NC_"
        nc_dict = {base_accver: (chrstart, chrstop) for base_accver, (chrstart, chrstop) in nc_dict.items() if
                   base_accver.startswith("NC_")}

        # Si le dictionnaire n'est pas vide, récupérez la base avant le point du premier élément
        if nc_dict:
            first_base = next(iter(nc_dict)).split('.')[0]  # Récupérer la base avant le point de la première clé

            # Filtrer le dictionnaire pour ne garder que les éléments avec la même base
            nc_dict = {base_accver: (chrstart, chrstop) for base_accver, (chrstart, chrstop) in nc_dict.items() if
                       base_accver.split('.')[0] == first_base}

        # Votre code existant pour afficher le dictionnaire filtré
        print("\nDictionnaire filtré :")
        for base_accver, (chrstart, chrstop) in nc_dict.items():
            print(f"{base_accver}: chrstart={chrstart}, chrstop={chrstop}")

        # Trouver l'entrée avec le chiffre le plus grand et le plus petit après le point
        max_version = -1
        max_accver = None
        max_coords = None  # Pour stocker les coordonnées associées à la version maximale
        min_version = float('inf')  # Initialiser à l'infini pour trouver la plus petite version
        min_accver = None
        min_coords = None  # Pour stocker les coordonnées associées à la version minimale

        for base_accver in nc_dict.keys():
            version = int(base_accver.split('.')[1])  # Extraire la version après le point

            # Mise à jour pour la version maximale
            if version > max_version:
                max_version = version
                max_accver = base_accver
                max_coords = nc_dict[base_accver]  # Stocker les coordonnées

            # Mise à jour pour la version minimale
            if version < min_version:
                min_version = version
                min_accver = base_accver
                min_coords = nc_dict[base_accver]  # Stocker les coordonnées

        # Afficher les résultats
        if max_accver:
            print(f"\nL'accès avec la version la plus élevée est : {max_accver} avec la version {max_version}.")
            print(f"Coordonnées : chrstart={max_coords[0]}, chrstop={max_coords[1]}")
        else:
            print("\nAucun accès trouvé.")

        if min_accver:
            print(f"L'accès avec la version la plus basse est : {min_accver} avec la version {min_version}.")
            print(f"Coordonnées : chrstart={min_coords[0]}, chrstop={min_coords[1]}")
        else:
            print("Aucun accès trouvé.")

        # Étape 3 : Faire une requête sur chaque accession version pour récupérer les infos détaillées
        for base_accver in nc_dict.keys():
            nc_info = fetch_nc_info(base_accver)
            if nc_info and 'result' in nc_info and base_accver in nc_info['result']:
                nc_details = nc_info['result'][base_accver]
                print(f"\nDétails pour {base_accver} :")
                print(f"Organisme : {nc_details.get('organism', {}).get('scientificname', 'Inconnu')}")
                print(f"Description : {nc_details.get('title', 'Inconnu')}")
            else:
                print(f"Aucune information détaillée trouvée pour {base_accver}.")

    else:
        print(f"Aucune information détaillée trouvée pour l'ID {gene_id}.")


if __name__ == "__main__":
    main()
