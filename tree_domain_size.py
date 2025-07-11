import csv


def read_treelist(filepath, newTreeFile, nx, ny):
    trees_list = []
    with open(filepath, "r") as file:
        reader = csv.reader(file, delimiter="\t")
        lines = list(reader)
        for line in lines:
            x_value = float(line[1])
            y_value = float(line[2])
            c_radius = float(line[4])
            if x_value + c_radius+1 >= nx or x_value-(c_radius+1) <= 0:
                print(f'x {x_value} with {c_radius} exceeds domain')
                continue
            if y_value+c_radius+1 >= nx or y_value-(c_radius+1) <= 0:
                print(f'y {y_value} with {c_radius} exceeds domain')
                continue
            trees_list.append(line)
    with open(newTreeFile, "w") as file:
        writer = csv.writer(file, delimiter="\t")
        for line in trees_list:
            writer.writerow(line)


if __name__ == "__main__":
    read_treelist("control/controlFireTec.txt", "control/controlClean.txt", 600, 600)


    
