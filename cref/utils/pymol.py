from pymol import cmd

cmd.show_as("cartoon", "experimental_structure")
cmd.show_as("cartoon", "predicted_structure")
cmd.do('''super predicted_structure, experimental_structure''')
