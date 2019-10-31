def read_target_plates(y):
        """ Reads the coordinates for a line defining the inner
        and outer target plates.
        The lines to define target plates end up with the shape ((x,y),(x,y)).
        These files read can contain more complicated lines than this.
        """
        target_plates = y['target_plates']
        d = {}
        for plate in target_plates:

            d[plate] = []
            if not 'name' in target_plates[plate].keys():
                target_plates[plate].update({'name' : plate})
                
            try:
                with open(target_plates[plate]['file']) as f:
                    for line in f:
                        point = line.strip()
                        if point.startswith('#'):
                            # move along to the next iteration
                            # this simulates a comment
                            continue
                        x = float(point.split(',')[0])
                        y = float(point.split(',')[1])
                        d[plate].append((x, y))
                print('Using target plate "{}": {}'.format(target_plates[plate]['name'], d[plate]))
            except KeyError:
                print('No inner target plate file.')


        return target_plates, d

def plot_target_plate(plate_data):
    import matplotlib
    matplotlib.use('TKagg')
    import matplotlib.pyplot as plt
    import numpy as np
    """ Plots the inner and outer target plates on the current figure """
    
    for plate in plate_data:
        coor = np.array(plate_data[plate])
        plt.plot(coor[:, 0], coor[:, 1], label=plate)
        plt.draw()

import yaml
import pathlib as p


f = p.Path('Plates_yaml')
y = yaml.load(f.read_text())

print(y['target_plates'].keys())

y['target_plates'], d = read_target_plates(y)

f.write_text(yaml.dump(y))

print('HERE IS DDDDDDDDDDDD: {}'.format(d))

for item in d:
    print(d[item])
    print(type(d[item]))

plot_target_plate(d)
