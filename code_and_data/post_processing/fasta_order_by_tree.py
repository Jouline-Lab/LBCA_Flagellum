import argparse
from ete3 import Tree

def order_fasta_based_on_tree(newick_tree, fasta_file):
    # Parse the Newick tree
    tree = Tree(newick_tree, format=1)  # format=1 assumes that leaf nodes are named

    # Extract the order of leaf nodes (sequence identifiers)
    leaf_order = [leaf.name for leaf in tree.iter_leaves()]

    # Read the FASTA file and store sequences in a dict keyed by their identifier
    fasta_dict = {}
    with open(fasta_file, 'r') as fasta:
        identifier = ""
        sequence = ""
        for line in fasta:
            line = line.strip()
            if line.startswith(">"):  # Header line
                if identifier:  # Save the previous sequence
                    fasta_dict[identifier] = sequence
                identifier = line[1:]  # Remove the '>' and store the new identifier
                sequence = ""  # Reset sequence for the new record
            else:
                sequence += line
        # Don't forget the last entry
        if identifier:
            fasta_dict[identifier] = sequence

    # Reorder the sequences based on the tree's leaf node order
    ordered_fasta = [(id, fasta_dict[id]) for id in leaf_order if id in fasta_dict]

    # Write the reordered sequences to a new file
    output_file = fasta_file.rsplit('.', 1)[0] + '_treeordered.fasta'  # Change file extension
    with open(output_file, 'w') as f:
        for id, seq in ordered_fasta:
            f.write(f">{id}\n{seq}\n")

    print(f"Reordered FASTA file written to {output_file}")

def main():
    parser = argparse.ArgumentParser(description="Order FASTA sequences based on a Newick tree")
    parser.add_argument('--tree', required=True, help="Path to the Newick tree file")
    parser.add_argument('--fasta', required=True, help="Path to the FASTA file")

    args = parser.parse_args()

    order_fasta_based_on_tree(args.tree, args.fasta)

if __name__ == "__main__":
    main()