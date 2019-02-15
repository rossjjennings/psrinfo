from subprocess import Popen, PIPE
import os.path

this_dir = os.path.dirname(__file__)
parent_dir = os.path.dirname(this_dir)
ne2001_binary_path = os.path.join(parent_dir, 'NE2001', 'bin.NE2001', 'NE2001')
ne2001_input_path = os.path.join(parent_dir, 'NE2001', 'input.NE2001')
ymw16_binary_path = os.path.join(parent_dir, 'ymw16', 'ymw16')
ymw16_input_path = os.path.join(parent_dir, 'ymw16', '')

def ne2001_dm_from_dist(pulsar, dist=None):
    '''
    Estimate the DM of `pulsar`, given its distance,
    using the NE2001 electron density model.
    
    Input:
    ------
    `pulsar`: Pulsar for which to estimate distance.
    `dist`: Distance, in kpc.
    
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
    args = [ne2001_binary_path, str(l), str(b), str(dist), '-1']
    process = Popen(args, stdout=PIPE, cwd=ne2001_input_path)
    stdout, stderr = process.communicate()
    stdout = stdout.decode('utf-8')
    words = stdout.split()
    return float(words[words.index('DM') - 1])

def ne2001_dist_from_dm(pulsar, dm=None):
    '''
    Estimate the distance to `pulsar` from its dispersion measure,
    using the NE2001 electron density model.
    
    Input:
    ------
    `pulsar`: Pulsar for which to estimate distance.
    `dm`: Dispersion measure, in pc/cm**3.
    
    Output:
    -------
    `dist`: Distance, in kpc.
    '''
    if dm is None:
        dm = pulsar.dm
    galactic_coords = pulsar.galactic()
    l = galactic_coords.l.value
    b = galactic_coords.b.value
    args = [ne2001_binary_path, str(l), str(b), str(dm), '1']
    process = Popen(args, stdout=PIPE, cwd=ne2001_input_path)
    stdout, stderr = process.communicate()
    stdout = stdout.decode('utf-8')
    words = stdout.split()
    return float(words[words.index('DIST') - 1])

def ymw16_dm_from_dist(pulsar, dist=None):
    '''
    Estimate the DM of `pulsar`, given its distance,
    using the YMW16 electron density model.
    
    Input:
    ------
    `pulsar`: Pulsar for which to estimate distance.
    `dist`: Distance, in kpc.
    
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
    args = [ymw16_binary_path, '-d', ymw16_input_path, 'Gal',
            str(l), str(b), str(dist*1000), '2']
    process = Popen(args, stdout=PIPE)
    stdout, stderr = process.communicate()
    stdout = stdout.decode('utf-8')
    words = stdout.split()
    return float(words[words.index('DM:') + 1])

def ymw16_dist_from_dm(pulsar, dm=None):
    '''
    Estimate the distance to `pulsar` from its dispersion measure,
    using the YMW16 electron density model.
    
    Input:
    ------
    `pulsar`: Pulsar for which to estimate distance.
    `dm`: Dispersion measure, in pc/cm**3.
    
    Output:
    -------
    `dist`: Distance, in kpc.
    '''
    if dm is None:
        dm = pulsar.dm
    galactic_coords = pulsar.galactic()
    l = galactic_coords.l.value
    b = galactic_coords.b.value
    args = [ymw16_binary_path, '-d', ymw16_input_path, 'Gal',
            str(l), str(b), str(dm), '1']
    process = Popen(args, stdout=PIPE)
    stdout, stderr = process.communicate()
    stdout = stdout.decode('utf-8')
    words = stdout.split()
    return float(words[words.index('Dist:') + 1])/1000
