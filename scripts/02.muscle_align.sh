muscle -align combined_monomer.fasta --output selected_combined_representative_monomer_align.fa
FastTreeMP -nt selected_combined_representative_monomer_align.fa > selected_combined_representative_monomer_align.tree

#muscle -align Zpal_monomer.fasta --output selected_Zpal_representative_monomer_align.fa
#FastTreeMP -nt selected_Zpal_representative_monomer_align.fa > selected_Zpal_representative_monomer_align.tree

#muscle -align Nip_monomer.fasta --output selected_Nip_representative_monomer_align.fa
#FastTreeMP -nt selected_Nip_representative_monomer_align.fa > selected_Nip_representative_monomer_align.tree
