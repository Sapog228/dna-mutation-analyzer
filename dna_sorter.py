codon_table = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}


def analyze_mutations(ref_dna: str, sample_dna: str, start_pos: int):
    """
    ref_dna, sample_dna - входные последовательности ДНК (референс и пример)
    start_pos -  позиция старта
    """

    # 1. Поиск мутаций (группируем подряд идущие пробелы в единые INS (вставки) / DEL (разделения))
    i = 0
    raw_variants = []  # (позиция мутации(pos), тип(type), символ/ы из последовательности референса(ref_seq), символ/ы из последовательности примера(alt_seq))

    if len(ref_dna) != len(sample_dna):
        raise ValueError("Ошибка: длины выровненных последовательностей должны совпадать.")

    n = len(ref_dna)

    if start_pos < 1 or start_pos > n:
        raise ValueError(f"Ошибка: позиция старта ({start_pos}) выходит за пределы длины выравнивания (1 - {n}).")

    valid_chars = set("ATGC-")
    if not set(ref_dna).issubset(valid_chars) or not set(sample_dna).issubset(valid_chars):
        raise ValueError("Ошибка: последовательности содержат недопустимые символы. Разрешены только A, T, G, C и '-'.")

    mutations = []
    frameshift_info = []

    while i < n:  # Берём попарно символы последовательностей ДНК (r - символ из ДНК референса, s - символ из ДНК примера)
        pos = i + 1  # Т.к. позиция должна быть не в виде индекса python(0-based), а в стандартном счете(1-based), то добавляем +1
        r, s = ref_dna[i], sample_dna[i]

        if r == s:  # Если символы равны, то идём дальше
            i += 1
            continue

        if r == '-':  # Если символы не равны и символы референса (r) это пробел, то
            ins_seq = []
            while i < n and ref_dna[i] == '-':
                ins_seq.append(sample_dna[
                                   i])  # Добавляем в лист символы из примера, которые стоят на местах пробела в референсе и объединяем их
                i += 1
            raw_variants.append((pos, 'Вставка', '-', "".join(
                ins_seq)))  # Добавляем в raw_variants данные, где произошла мутация, её тип, на месте какого символа она произошла, на какие символы заменилась

        elif s == '-':  # Если символы не равны и символ примера (s) это пробел, то
            del_seq = []
            while i < n and sample_dna[i] == '-':
                del_seq.append(ref_dna[
                                   i])  # Добавляем в лист символ/ы, которые стояли на местах пробелов в референсе и объединяем их
                i += 1
            raw_variants.append((pos, 'Делеция', "".join(del_seq),
                                 '-'))  # Добавляем в raw_variants данные, где произошла мутация, её тип, на месте каких символов она произошла, на какие символы заменилась

        else:  # Одиночная замена (SNP), буквы просто разные
            raw_variants.append((pos, 'Замена', r,
                                 s))  # Добавляем в raw_variants данные, где произошла мутация, её тип, на месте какого символов она произошла, на какой символ заменилась
            i += 1

    # Вспомогательные карты индексов (от старта рамки считывания) для точного маппинга без '-'
    ref_coding_indices = [idx for idx in range(start_pos - 1, n) if ref_dna[idx] != '-']
    sample_coding_indices = [idx for idx in range(start_pos - 1, n) if sample_dna[idx] != '-']

    # 2. Оценка эффекта каждой мутации
    for pos, var_type, ref_al, alt_al in raw_variants:  # считываем позицию(pos), тип мутации(var_type), правильный аллель(строку) из референса(ref_al), алелль с мутацией(alt_al)
        var_len = len(alt_al) if var_type == 'Вставка' else len(ref_al)
        var_end_pos = pos + var_len - 1

        # Проверка на некодирующую область (до start_pos), если мутация произошла до стартовой точки кодирования ДНК
        if var_end_pos < start_pos:
            mutations.append({
                'Позиция': pos,
                'Тип': var_type,
                'Реф. аллель': ref_al,
                'Прим. аллель': alt_al,
                'Эффект': 'Некодирующий'  # Означает, что не влияет на ДНК напрямую, но обязательно к записи
            })
            continue

        if var_type == 'Замена':
            # Смещение внутри кодирующей рамки (0-based)
            aln_idx = pos - 1  # Смещение, на сколько букв от старта стоит мутация
            ref_nt_offset = ref_coding_indices.index(aln_idx)
            codon_idx = ref_nt_offset // 3  # (offset // 3) необходимо, чтобы найти в каком кодоне(строке из 3х символов) находится мутация
            codon_pos_in_triplet = ref_nt_offset % 3  # Нахождение нужного символа в кодоне

            codon_start_nt = codon_idx * 3
            if codon_start_nt + 3 > len(ref_coding_indices):
                effect = 'Некодирующий'
            else:
                ref_codon_idxs = ref_coding_indices[codon_start_nt:codon_start_nt + 3]
                ref_codon = "".join(ref_dna[idx] for idx in ref_codon_idxs)  # Копирование кодона из ДНК референса

                alt_chars = [sample_dna[idx] for idx in ref_codon_idxs]
                alt_chars[codon_pos_in_triplet] = alt_al
                alt_codon = "".join(alt_chars)  # Копирование кодона из ДНК примера

                if '-' in alt_codon:  # Если в кодон попал пробел, то невозможно его декодировать
                    effect = 'Сдвиг рамки считывания'
                else:
                    ref_aa = codon_table.get(ref_codon)  # Сравнение кодона с переводом из словаря
                    alt_aa = codon_table.get(alt_codon)  # Сравнение кодона с переводом из словаря

                    if ref_aa is None or alt_aa is None:  # Если не возможно перевести хотя бы один кодон, то мы эффект не известен
                        effect = 'Неизвестно'
                    elif ref_aa == alt_aa:
                        effect = 'Синонимичный'  # Если перевод кодонов идентичный
                    else:
                        effect = 'Несинонимичный'  # Если перевод кодонов различается

            mutations.append({
                'Позиция': pos,
                'Тип': var_type,
                'Реф. аллель': ref_al,
                'Прим. аллель': alt_al,
                'Эффект': effect
            })

        else:  # Когда тип мутации вставка(INS) или делеция(DEL)
            indel_len = len(alt_al) if var_type == 'Вставка' else len(
                ref_al)  # InDel (INS + DEL). Берётся длина аллели примера при INS и длина реф. аллели при DEL

            if indel_len % 3 != 0:  # Длина не делится цело на 3, то это эффект сдвига рамки считывания
                effect = 'Сдвиг рамки считывания'

                # Определение индекса первого затронутого кодона (0-based от старта трансляции)
                # Расчёт, кол-ва нуклеотидов(символа перевода кодона) в sample до позиции мутации
                sample_nt_before = sum(1 for idx in sample_coding_indices if idx < pos - 1)
                affected_codon_idx = sample_nt_before // 3  # Расчёт кодона, где произошла мутация

                total_sample_codons = len(sample_coding_indices) // 3

                # Транслируем мутантную рамку от точки сбоя до стоп-кодона
                affected_aa_count = 0
                new_stop_pos = None  # Т.к. из-за мутации может образоваться стоп-кодон(например: TGA: *), то сборка ДНК может закончиться раньше времени и необходимо это учесть.

                for c_idx in range(affected_codon_idx, total_sample_codons):  # Расчёт начинается с проблемного кодона
                    nt_slice = sample_coding_indices[c_idx * 3: c_idx * 3 + 3]
                    codon = "".join(sample_dna[idx] for idx in nt_slice)  # Срез кодона.
                    aa = codon_table.get(codon, '?')  # Перевод кодона в аминокислоту

                    if aa == '*':
                        # Координата стоп-кодона в исходной кодирующей последовательности (нуклеотидная)
                        new_stop_pos = nt_slice[0] + 1
                        break
                    affected_aa_count += 1

                frameshift_info.append((pos, new_stop_pos,
                                        affected_aa_count))  # (Позиция(pos), новая_стоп_позиция(если есть(new_stop_pos)), количество_затронутых_аминокислот(affected_aa_count))

            else:
                # Длина кратна 3 — сдвига рамки нет
                effect = 'Без сдвига рамки'

            mutations.append({
                'Позиция': pos,
                'Тип': var_type,
                'Реф. аллель': ref_al,
                'Прим. аллель': alt_al,
                'Эффект': effect
            })

    return mutations, frameshift_info


if __name__ == "__main__":
    ref = "GCATCGATCG--CGATCGATTTATCGATCGATGGCAAGTATTCCTAAGTTAAAGGGCCCTATTGCGATGAACACCAACGTATTACTAATGTTGCT---GTGACAGCTCACAGCTGGTCGCATCGATCACCACAACAAGTCGCCGGCGTTTTTATCCCTACTGCCTCCCCCACCGCATCACCAACGTCGCCGGCGATGA---TCTACCACAACGGATAACAAAAGTAGCAGGAGACGAGGGGCTGCCGCAGCGGCGGATGACGAAGGTGGCGGGGGGATGAAGGTCTACCACAACGAATAACAAAAGTAGCAGGAGACGAGGGGCTGCCGCAGCGGATGACGAAGGTGGCGGGGGATGAAGGTCTACCACAACGAATAACAAAGTAGCAGGAGACGAGGGGCTGCCGCAGCGGATGACGAAGGTGGCGGGGGATGAAGGTCTACCACAACGAATAACAAAGTAGCAGGAGACGAGGGGCTGCCGTAATTTGGGCCCAAAATTGCCCCCTTTTA"
    sample = "GCATAGATCGATCGATCGAT--ATCGATCGATGGCAAGAATTC---AGTTAAAGGGCCCTATTGCGATGAACACCAACGTATTACTAATGTTGCTGCTGTGACAGCTCACAGCTG-CCGCATCGATCACCACAACAAGTCGCCGGCGTTTTTATCCCTAC---CTCCCCCACCGCATCACCAACGTCGCCGGCGATGAAGGTCTACCACAACGGATAACAAAAGTAGCAGGAGACGAGGGGCTGCCGCAGCGGCGGATGACGAAGGTGGCGGGGGGATGAAGGTCTACCACAACGAATAACAAAAGTAGCAGGAGACGAGGAAAAACCGCAGCGGATGACGAAGGTGGCGGGGGATGAAGGTCTACCACAACGAATAACAAAGTAGCAGGAGACGAGGGGCTGCCGCAGCGGATGACGAAGGTGGCGGGGGATGAAGGTCTACCACAACGAATAACAAAGTAGCAGGAGACGAGGGGCTGCCGTAATTTGGGCCCAATTTTGCCCCCTTTTA"
    start = 31

    muts, fs_info = analyze_mutations(ref, sample, start)

    print("=== ТАБЛИЦА МУТАЦИЙ ===")
    header = f"{'Позиция':<10} | {'Тип':<10} | {'Реф':<6} | {'Прим':<6} | {'Эффект'}"
    print(header)
    print("-" * len(header))

    for m in muts:
        print(f"{m['Позиция']:<10} | {m['Тип']:<10} | {m['Реф. аллель']:<6} | {m['Прим. аллель']:<6} | {m['Эффект']}")

    print("\n=== СТАТИСТИКА СДВИГА РАМКИ ===")
    if not fs_info:
        print("Мутаций со сдвигом рамки нет.")
    else:
        fs_header = f"{'Позиция мутации':<18} | {'Новый стоп-кодон':<18} | {'Затронуто аминокислот'}"
        print(fs_header)
        print("-" * len(fs_header))
        for p, stop_p, count in fs_info:
            stop_display = str(stop_p) if stop_p is not None else "Не найден"
            print(f"{p:<18} | {stop_display:<18} | {count}")