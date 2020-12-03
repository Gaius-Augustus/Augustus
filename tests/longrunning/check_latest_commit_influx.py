#!/usr/bin/env python3

import sys
import argparse
from influxdb import InfluxDBClient

# Output the hash of the last tested commit from the influxdb.

parser = argparse.ArgumentParser(
    description='Get latest commit stored in InfluxDB.')
parser.add_argument('-d', '--dbname', help='name of the database to use.')
parser.add_argument('-o', '--dbhost', help='host of influxdb.')
parser.add_argument('-p', '--dbport', help='port of influxdb.')
parser.add_argument('-u', '--dbuser', help='user of the database to use.')
parser.add_argument('-w', '--dbpasswd',
                    help='password of the given database user.')
args = parser.parse_args()


def check_last_commit():
    client = InfluxDBClient(host=args.dbhost, port=args.dbport,
                            username=args.dbuser, password=args.dbpasswd, database=args.dbname)

    res = client.query(
        "SELECT * FROM qualitymetrics WHERE mode='cgp' ORDER BY time DESC LIMIT 1")
    latest_hash = list(res.get_points(measurement='qualitymetrics'))[0]['hash']
    return latest_hash


if __name__ == '__main__':
    if args.dbname is None:
        print('The database to use is required, please make use of --dbname to pass the name of the database...')
        sys.exit()

    if args.dbhost is None:
        print('The host of influxdb is required, please make use of --dbhost to pass the host...')
        sys.exit()

    if args.dbport is None:
        print('The port of influxdb is required, please make use of --dbport to pass the name of the database...')
        sys.exit()

    if args.dbuser is None:
        print(
            'The database user is required, please make use of --dbuser to pass the user...')
        sys.exit()

    if args.dbpasswd is None:
        print('The password of the given user is required, please make use of --dbpasswd to pass the password...')
        sys.exit()

    print(check_last_commit())
