import csv
import collections

# Dwyer Data from doi:10.1186/1478-7954-12-5
dwyer = list()
with open('Dwyer-S3.txt') as read_file:
    reader = csv.DictReader(read_file, delimiter='\t')
    for row in reader:
        for year in range(1996, 2013) + ['Annualized rate of change, 1996-2012']:
            year = str(year)
            prevalence, prevalence_CI = row[year].split(' ')
            lower, upper = prevalence_CI.strip('()').split(',')
            row[year] = float(prevalence)
            row[year + '_lower'], row[year + '_upper'] = float(lower), float(upper)
        if not row['County']:
            row['County'] = dwyer[-1]['County']
        dwyer.append(row)

for row in dwyer:
    row['County'] = row['County'].replace(' County', '', 1)

dwyer = [row for row in dwyer if row['Sex'] == 'Both']

with open('lung-table.txt') as read_file:
    reader = csv.DictReader(read_file, delimiter='\t')
    lung_fields = reader.fieldnames
    lung = list(reader)

county_to_dwyer = {county['County'].lower(): county for county in dwyer}

for row in lung:
    name = row['name'].lower()
    county = county_to_dwyer.get(name)
    row['smoking_1996'] = county['1996'] if county else None
    row['smoking_2012'] = county['2012'] if county else None
    row['dwyer_delta'] = county['Annualized rate of change, 1996-2012'] if county else None

print collections.Counter(bool(row['smoking_1996']) for row in lung)

with open('lung-table-with-dwyer.txt', 'w') as write_file:
    fieldnames = lung_fields + ['smoking_1996', 'smoking_2012', 'dwyer_delta']
    writer = csv.DictWriter(write_file, delimiter = '\t', fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(lung)

