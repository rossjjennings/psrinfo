import numpy as np
from subprocess import Popen, PIPE
from io import BytesIO
from collections import namedtuple
import os.path
from ruamel.yaml import YAML
yaml = YAML(typ='safe')

from psrinfo.psr import pulsar

this_dir = os.path.dirname(__file__)
parent_dir = os.path.dirname(this_dir)
psrcat_binary_path = os.path.join(parent_dir, 'psrcat', 'psrcat')
psrcat_db_path = os.path.join(parent_dir, 'psrcat', 'psrcat.db')

with open(os.path.join(this_dir, 'psrcat_params.yml')) as infile:
    psrcat_params = yaml.load(infile)

def fetch_pulsar(name, extra_params=None):
    record = fetch_records(extra_args=[name], extra_params=extra_params)
    pulsar_tuple = namedtuple('pulsar_tuple', record.dtype.names)
    t = pulsar_tuple(*record.item())
    if len(t.elat_cite) > 1:
        return pulsar.from_ecliptic(t)
    else:
        return pulsar.from_equatorial(t)

def fetch_pulsars(condition=None, extra_params=None):
    extra_args = list()
    if condition is not None:
        extra_args += ['-l', condition]
    records = fetch_records(extra_args=extra_args, extra_params=extra_params)
    pulsar_tuple = namedtuple('pulsar_tuple', records.dtype.names)
    tuples = [pulsar_tuple(*row) for row in records]
    pulsars = {t.name: pulsar.from_ecliptic(t) 
               if len(t.elat_cite) > 1 else pulsar.from_equatorial(t)
               for t in tuples}
    return pulsars

def fetch_records(extra_args=None, extra_params=None):
    args = [psrcat_binary_path, '-db_file', psrcat_db_path]
    args += ['-o', 'long_error_csv']
    if extra_args is not None:
        args += extra_args
    params = ['name', 'raj', 'decj', 'elat', 'elong', 'dm',
              'pmra', 'pmdec', 'pmelat', 'pmelong', 'posepoch']
    if extra_params is not None:
        for param in extra_params:
            params.append(param)
    args += ['-c', ' '.join(params)]
    
    fields = generate_fields(params)
    process = Popen(args, stdout=PIPE)
    stdout, stderr = process.communicate()
    records = np.recfromtxt(BytesIO(stdout), encoding='utf-8',
                            skip_header=2, delimiter=';',
                            usecols=range(len(fields)), names=fields)
    return records

def generate_fields(params):
    fields = ['num']
    for param in params:
        if param.upper() not in psrcat_params:
            error_msg = "'{}' is not a recognized parameter name"
            raise ValueError(error_msg.format(param))
        if param.upper() in ['RAJ', 'DECJ', 'ELONG', 'PMELONG']:
            fieldname = param[:-1].lower()
        else:
            fieldname = param.lower()
        fields.append(fieldname)
        if psrcat_params[param.upper()] in ['dec', 'sexg']:
            fields.append(fieldname + '_err')
        if param.upper() == 'NAME':
            fields.append('cite')
        else:
            fields.append(fieldname + '_cite')
    return fields
