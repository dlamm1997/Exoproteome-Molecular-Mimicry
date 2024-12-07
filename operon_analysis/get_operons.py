import csv
import sys



def find_operon_by_gene_id(file_name, gene_id):
    try:
        with open(file_name, 'r', newline='', encoding='utf-8') as file:
            reader = csv.DictReader(file, delimiter='\t')
	    print(reader)
            for row in reader:
                if row['IdGene'] == gene_id:
                    return row['Operon']
        return f"Gene ID {gene_id} not found in the file."
    except FileNotFoundError:
        return f"Error: File {file_name} not found."
    except KeyError:
        return "Error: One of the expected columns is missing in the file."
    except Exception as e:
        return str(e)

if __name__ == "__main__":
    # Example usage
    file_name_input = sys.argv[2]
    gene_id_input = sys.argv[1]
    
    operon_result = find_operon_by_gene_id(file_name_input, gene_id_input)
    print(f"The operon associated with IdGene {gene_id_input} is: {operon_result}")
