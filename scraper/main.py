import sys
from queue import Queue

sys.path.append('..')
from scraper.get import *


def main():
    # input
    root_url = 'http://ilthermo.boulder.nist.gov'
    init_ion = 'fluoride'
    # First we get property table
    prp_table = get_prp_table(root_url + '/ILT2/ilprpls')
    prp_index = put_prp_table(prp_table)

    search_queue = Queue()
    exist_ions = [ion.name for ion in session.query(Ion).filter_by(searched=False)]
    if exist_ions:
        list(map(search_queue.put, exist_ions))
    else:
        search_queue.put(init_ion)
        ion = put_unique_ion(init_ion)
        put_ion(init_ion, ion)

    session.commit()
    idxcut = -1
    while not search_queue.empty():
        print('There are %d species in queue' % search_queue.qsize())
        search_name = search_queue.get()
        print('Search species:', search_name)

        try:
            paper_table = get_paper_table(root_url + '/ILT2/ilsearch', params={
                'cmp': search_name,
                'ncmp': '',
                'year': '',
                'auth': '',
                'keyw': '',
                'prp': ''
            })
        except SearchFailedError:
            Log.write('Cannot search. Skip this...')
            continue
        except SpecialCaseError:
            Log.write('Special case error. Skip this...')
            continue

        print('Get %d papers, ' % len(paper_table), end='')
        paper_table_new = []
        for line in paper_table:
            if not session.query(DataSet).filter_by(code=line['code']).first():
                paper_table_new.append(line)
        paper_table = paper_table_new
        print('%d need to be downloaded' % len(paper_table))

        for idx, line in enumerate(paper_table):
            print('idx=', idx, line['code'])
            if line['code'] in ['pxEeK', 'BExHp', 'rjzQF', 'leTaX', 'lfBxY', 'PwbXr']:
                continue
            # print('[%d%%] Search paper %s %s (%s)...' % (idx*100/len(paper_table), line['code'], line['property'],
            # line['ref']), end='')
            if idx < idxcut:
                continue
            try:
                paper_info, molecule_info, raw_data = get_paper_molecule_info(root_url + '/ILT2/ilset',
                                                                              params={'set': line['code']})
            except SearchFailedError:
                Log.write('Cannot get data from paper. Skipping...')
                continue
            except SpecialCaseError:
                Log.write('Cannot read data from paper properly. Skipping...')
                continue
            put_paper(paper_info)
            for m_info in molecule_info:
                if session.query(Molecule).filter_by(name=m_info.get('name')).first():
                    print('Molecule:', m_info.get('name'), 'exists in database, id=',
                          session.query(Molecule).filter_by(name=m_info.get('name')).first().id)
                    continue
                else:
                    print('Molecule:', m_info.get('name'))
                cation = anion = None
                if m_info.get('cation') is not None:
                    cation_name = m_info.get('cation')
                    cation = put_unique_ion(cation_name)
                    if cation.charge <= 0:
                        raise Exception('cation charge error')
                    if not session.query(Ion).filter_by(name=cation_name).first():
                        search_queue.put(cation_name)
                        put_ion(cation_name, cation)
                if m_info.get('anion') is not None:
                    anion_name = m_info.get('anion')
                    anion = put_unique_ion(anion_name)
                    if anion.charge >= 0:
                        raise Exception('anion charge error')
                    if not session.query(Ion).filter_by(name=anion_name).first():
                        search_queue.put(anion_name)
                        put_ion(anion_name, anion)
                # print(cation, anion, type(cation), type(anion))
                mol = put_unique_molecule(m_info, cation, anion)
                put_molecule(m_info, mol)
            session.add(DataSet(code=line['code'], raw_data=raw_data))
            session.commit()
            print('')
        session.query(Ion).filter_by(name=search_name).update({Ion.searched: True})
        session.commit()


if __name__ == '__main__':
    main()
