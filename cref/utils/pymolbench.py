import tempfile
from pymol import cmd

cmd.show_as("cartoon", "experimental_structure")
cmd.show_as("cartoon", "predicted_structure")
rmsd = cmd.align('predicted_structure', 'experimental_structure')
cmd.bg_color('white')
handler, output_file = tempfile.mkstemp(prefix='alignment', suffix='.png')
cmd.png(output_file, ray=1)
print(rmsd[0])
print(output_file)
cmd.quit()
