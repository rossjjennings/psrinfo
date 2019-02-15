from subprocess import Popen, PIPE
import os.path

this_dir = os.path.dirname(__file__)
parent_dir = os.path.dirname(this_dir)
ne2001_binary_path = os.path.join(parent_dir, 'NE2001', 'bin.NE2001', 'NE2001')
ne2001_input_path = os.path.join(parent_dir, 'NE2001', 'input.NE2001')
ymw16_binary_path = os.path.join(parent_dir, 'ymw16', 'ymw16')
ymw16_input_path = os.path.join(parent_dir, 'ymw16', '')

def dist_from_dm(pulsar, model='NE2001', dm=None):
    '''
    Estimate the DM of `pulsar`, given its distance,
    using a galactic electron density model.
    
    Input:
    ------
    `pulsar`: Pulsar for which to estimate distance.
    `dm`: Dispersion measure, in pc/cm**3.
    `model`: Electron density model to use, one of 'NE2001' or 'YMW16'.
    
    Output:
    -------
    `dist`: Distance, in kpc.
    '''
    if dm is None:
        dm = pulsar.dm
    galactic_coords = pulsar.galactic()
    l = galactic_coords.l.value
    b = galactic_coords.b.value
    if model.upper() == 'NE2001':
        args = [ne2001_binary_path, str(l), str(b), str(dm), '1']
        process = Popen(args, stdout=PIPE, cwd=ne2001_input_path)
        stdout, stderr = process.communicate()
        stdout = stdout.decode('utf-8')
        words = stdout.split()
        dist = float(words[words.index('DIST') - 1])
    elif model.upper() == 'YMW16':
        args = [ymw16_binary_path, '-d', ymw16_input_path, 'Gal',
            str(l), str(b), str(dm), '1']
        process = Popen(args, stdout=PIPE)
        stdout, stderr = process.communicate()
        stdout = stdout.decode('utf-8')
        words = stdout.split()
        dist = float(words[words.index('Dist:') + 1])/1000
    return dist

def dm_from_dist(pulsar, model='NE2001', dist=None):
    '''
    Estimate the DM of `pulsar`, given its distance,
    using a galactic electron density model.
    
    Input:
    ------
    `pulsar`: Pulsar for which to estimate distance.
    `dist`: Distance, in kpc.
    `model`: Electron density model to use, one of 'NE2001' or 'YMW16'.
    
    Output:
    -------
    `dm`: Dispersion measure, in pc/cm**3.
    '''
    if dist is None:
        px = pulsar.px
        dist = 1/px
    galactic_coords = pulsar.galactic()
    l = galactic_coords.l.value
    b = galactic_coords.b.value
    if model.upper() == 'NE2001':
        args = [ne2001_binary_path, str(l), str(b), str(dist), '-1']
        process = Popen(args, stdout=PIPE, cwd=ne2001_input_path)
        stdout, stderr = process.communicate()
        stdout = stdout.decode('utf-8')
        words = stdout.split()
        dm = float(words[words.index('DM') - 1])
    elif model.upper() == 'YMW16':
        args = [ymw16_binary_path, '-d', ymw16_input_path, 'Gal',
            str(l), str(b), str(dist*1000), '2']
        process = Popen(args, stdout=PIPE)
        stdout, stderr = process.communicate()
        stdout = stdout.decode('utf-8')
        words = stdout.split()
        return float(words[words.index('DM:') + 1])
