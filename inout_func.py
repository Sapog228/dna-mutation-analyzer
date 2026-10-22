import csv
from collections import Counter
from dna_sorter import analyze_mutations

def load_from_csv(filepath):
    """Читает исходный CSV со штаммами"""
    data = []
    for enc in ('utf-8-sig', 'cp1251'):
        try:
            with open(filepath, newline='', encoding=enc) as f:
                reader = csv.DictReader(f, delimiter=';')
                for row in reader:
                    raw_name = row['strain_name'].strip()
                    data.append({
                        'strain_name': raw_name,
                        'ref_dna': row['ref_dna'].strip().upper(),
                        'sample_dna': row['sample_dna'].strip().upper(),
                        'start_pos': int(row['start_pos'])
                    })
            break  # Если успешно прочитали — выходим из цикла подбора кодировки
        except UnicodeDecodeError:
            data.clear()
            continue

    return data

def save_mutations_to_csv(all_records, output_filepath="mutations_report.csv"):
    """Записывает итоговый отчет обо всех мутациях в один CSV"""
    fieldnames = [
        'Штамм',
        'Эталон ДНК',
        'ДНК образца',
        'Позиция',
        'Тип',
        'Реф. аллель',
        'Прим. аллель',
        'Эффект',
        'Новый стоп-кодон',
        'Затронуто аминокислот'
    ]

    with open(output_filepath, mode='w', newline='', encoding='utf-8-sig') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter=';')
        writer.writeheader()
        writer.writerows(all_records)

    print(f"Таблица сохранена в {output_filepath}")


def print_statistics(records):
    """Считает статистику по мутациям"""
    print(f"\nСТАТИСТИКА МУТАЦИЙ:")
    print(f"Всего мутаций: {len(records)}")

    types_count = Counter(r['Тип'] for r in records)
    print("\nПо типу:")
    for t, c in types_count.items():
        print(f"  {t}: {c}")

    effects_count = Counter(r['Эффект'] for r in records)
    print("\nПо влиянию:")
    for e, c in effects_count.items():
        print(f"  {e}: {c}")

    fs_strains = {
        r['Штамм']
        for r in records
        if r['Эффект'] == 'Сдвиг рамки считывания'
    }
    print("\nШтаммы со сдвигом рамки считывания:")
    if not fs_strains:
        print("  Нет таких")
    else:
        for name in fs_strains:
            print(f"  {name}")


if __name__ == "__main__":
    input_csv = "strains_data.csv"
    output_csv = "mutations_report.csv"

    try:
        strains = load_from_csv(input_csv)
    except FileNotFoundError:
        print(f"Файл {input_csv} не найден")
        exit(1)

    all_mutation_rows = []

    for item in strains:
        name = item['strain_name']
        ref = item['ref_dna']
        sample = item['sample_dna']
        print(f"Обработка: {name}")

        muts, fs_info = analyze_mutations(
            ref,
            sample,
            item['start_pos']
        )

        fs_map = {pos: (stop_pos, count) for pos, stop_pos, count in fs_info}

        for m in muts:
            pos = m['Позиция']
            stop_pos, affected_aa = fs_map.get(pos, ('-', '-'))

            all_mutation_rows.append({
                'Штамм': name,
                'Эталон ДНК': ref,
                'ДНК образца': sample,
                'Позиция': pos,
                'Тип': m['Тип'],
                'Реф. аллель': m['Реф. аллель'],
                'Прим. аллель': m['Прим. аллель'],
                'Эффект': m['Эффект'],
                'Новый стоп-кодон': stop_pos if stop_pos is not None else 'Не найден',
                'Затронуто аминокислот': affected_aa
            })

    save_mutations_to_csv(all_mutation_rows, output_csv)
    print_statistics(all_mutation_rows)