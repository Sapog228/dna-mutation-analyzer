import tkinter as tk
from tkinter import ttk, filedialog, messagebox
from collections import Counter

from dna_sorter import analyze_mutations
from inout_func import load_from_csv, save_mutations_to_csv


class GradientHeader(tk.Canvas):
    """Шапка с градиентом в пастельно-голубых тонах"""

    def __init__(self, parent, color1="#4B86E8", color2="#8DC5F8", height=64, **kwargs):
        super().__init__(parent, height=height, highlightthickness=0, **kwargs)
        self.color1 = color1
        self.color2 = color2
        self.bind("<Configure>", self._draw_gradient)

    def _draw_gradient(self, event=None):
        self.delete("gradient")
        width = self.winfo_width()
        height = self.winfo_height()
        if width <= 1:
            return

        r1, g1, b1 = self.winfo_rgb(self.color1)
        r2, g2, b2 = self.winfo_rgb(self.color2)
        r1, g1, b1 = r1 // 256, g1 // 256, b1 // 256
        r2, g2, b2 = r2 // 256, g2 // 256, b2 // 256

        step = 2
        for i in range(0, width, step):
            factor = i / width
            nr = int(r1 + (r2 - r1) * factor)
            ng = int(g1 + (g2 - g1) * factor)
            nb = int(b1 + (b2 - b1) * factor)
            hex_color = f"#{nr:02x}{ng:02x}{nb:02x}"
            self.create_rectangle(i, 0, i + step, height, fill=hex_color, outline=hex_color, tags="gradient")

        self.tag_lower("gradient")


class DNAMutationApp(tk.Tk):
    COLOR_BG = "#F4F7FB"
    COLOR_CARD = "#FFFFFF"
    COLOR_PRIMARY = "#4A86E8"
    COLOR_PRIMARY_HOVER = "#3B73CE"
    COLOR_SEC = "#E6EFFC"
    COLOR_SEC_HOVER = "#D4E3FA"
    COLOR_TEXT = "#202E43"
    COLOR_MUTED = "#6F829D"
    COLOR_ROW_EVEN = "#FFFFFF"
    COLOR_ROW_ODD = "#EEF4FC"
    COLOR_BORDER = "#DCE4F0"

    def __init__(self):
        super().__init__()
        self.title("DNA Mutation Analyzer")
        self.geometry("1100x680")
        self.minsize(920, 520)
        self.configure(bg=self.COLOR_BG)

        self.input_file = ""
        self.all_results = []  # Все результаты целиком
        self.displayed_results = []  # Отфильтрованные данные

        self._init_styles()
        self._build_ui()

    def _init_styles(self):
        style = ttk.Style()
        style.theme_use("clam")

        # Таблица
        style.configure("Treeview",
                        background=self.COLOR_ROW_EVEN,
                        foreground=self.COLOR_TEXT,
                        fieldbackground=self.COLOR_ROW_EVEN,
                        rowheight=28,
                        font=("Segoe UI", 9),
                        borderwidth=0)

        style.configure("Treeview.Heading",
                        background="#E3EDF9",
                        foreground=self.COLOR_TEXT,
                        font=("Segoe UI", 9, "bold"),
                        relief="flat",
                        padding=6)

        style.map("Treeview",
                  background=[("selected", "#CCE2FD")],
                  foreground=[("selected", "#102A43")])

        style.map("Treeview.Heading",
                  background=[("active", "#D4E4F5")])

        # Кастомный современный скроллбар без громоздких 3D-рамок
        style.configure("Vertical.TScrollbar",
                        gripcount=0,
                        background="#C6DCF8",
                        troughcolor="#F4F7FB",
                        bordercolor="#F4F7FB",
                        darkcolor="#C6DCF8",
                        lightcolor="#C6DCF8",
                        arrowcolor="#4A86E8",
                        arrowsize=12,
                        relief="flat",
                        width=12)

        style.map("Vertical.TScrollbar",
                  background=[("active", "#A2C7F7"), ("pressed", "#7EAFF2")])

        # Выпадающие списки (Combobox)
        style.configure("Filter.TCombobox",
                        fieldbackground="#FFFFFF",
                        background="#E6EFFC",
                        foreground=self.COLOR_TEXT,
                        arrowcolor="#4A86E8",
                        relief="flat",
                        padding=4)
        style.map("Filter.TCombobox",
                  fieldbackground=[("readonly", "#FFFFFF")],
                  selectbackground=[("readonly", "#E6EFFC")],
                  selectforeground=[("readonly", self.COLOR_TEXT)])

    def _make_button(self, parent, text, cmd, primary=False, state=tk.NORMAL):
        bg = self.COLOR_PRIMARY if primary else self.COLOR_SEC
        fg = "#FFFFFF" if primary else self.COLOR_TEXT
        hover_bg = self.COLOR_PRIMARY_HOVER if primary else self.COLOR_SEC_HOVER

        btn = tk.Button(
            parent,
            text=text,
            command=cmd,
            bg=bg,
            fg=fg,
            activebackground=hover_bg,
            activeforeground=fg,
            font=("Segoe UI", 9, "bold" if primary else "normal"),
            bd=0,
            padx=14,
            pady=6,
            relief=tk.FLAT,
            cursor="hand2",
            state=state
        )

        def on_enter(e):
            if btn["state"] == tk.NORMAL:
                btn.configure(bg=hover_bg)

        def on_leave(e):
            if btn["state"] == tk.NORMAL:
                btn.configure(bg=bg)

        btn.bind("<Enter>", on_enter)
        btn.bind("<Leave>", on_leave)
        return btn

    def _build_ui(self):
        # 1. Шапка
        self.header = GradientHeader(self, color1="#4B86E8", color2="#8DC5F8", height=65)
        self.header.pack(fill=tk.X)

        self.header.create_text(25, 23, text="Анализатор Мутаций ДНК", anchor=tk.W,
                                fill="#FFFFFF", font=("Segoe UI", 15, "bold"))
        self.header.create_text(25, 45, text="Пакетный анализ последовательностей и фильтрация мутаций", anchor=tk.W,
                                fill="#D9EBFF", font=("Segoe UI", 9))

        # 2. Панель управления и фильтрации
        toolbar = tk.Frame(self, bg=self.COLOR_CARD, padx=15, pady=10,
                           highlightthickness=1, highlightbackground=self.COLOR_BORDER)
        toolbar.pack(fill=tk.X, padx=15, pady=(12, 6))

        # Левая часть — кнопки загрузки
        self.btn_select = self._make_button(toolbar, "Выбрать CSV", self.select_file)
        self.btn_select.pack(side=tk.LEFT)

        self.lbl_file = tk.Label(toolbar, text="Файл не выбран", bg=self.COLOR_CARD,
                                 fg=self.COLOR_MUTED, font=("Segoe UI", 9, "italic"))
        self.lbl_file.pack(side=tk.LEFT, padx=12)

        self.btn_run = self._make_button(toolbar, "Запустить анализ", self.run_analysis, primary=True,
                                         state=tk.DISABLED)
        self.btn_run.pack(side=tk.LEFT, padx=5)

        # Правая часть — экспорт
        self.btn_export = self._make_button(toolbar, "Экспорт отчета", self.export_csv, state=tk.DISABLED)
        self.btn_export.pack(side=tk.RIGHT)

        # Центр/Разделитель — фильтры
        filter_box = tk.Frame(toolbar, bg=self.COLOR_CARD)
        filter_box.pack(side=tk.RIGHT, padx=25)

        # Фильтр по типу
        tk.Label(filter_box, text="Тип:", bg=self.COLOR_CARD, fg=self.COLOR_MUTED,
                 font=("Segoe UI", 9, "bold")).pack(side=tk.LEFT, padx=(0, 4))
        self.type_var = tk.StringVar(value="Все")
        self.cb_type = ttk.Combobox(filter_box, textvariable=self.type_var, state="disabled",
                                    width=11, style="Filter.TCombobox")
        self.cb_type['values'] = ["Все"]
        self.cb_type.bind("<<ComboboxSelected>>", lambda e: self.apply_filters())
        self.cb_type.pack(side=tk.LEFT, padx=(0, 15))

        # Фильтр по эффекту
        tk.Label(filter_box, text="Эффект:", bg=self.COLOR_CARD, fg=self.COLOR_MUTED,
                 font=("Segoe UI", 9, "bold")).pack(side=tk.LEFT, padx=(0, 4))
        self.effect_var = tk.StringVar(value="Все")
        self.cb_effect = ttk.Combobox(filter_box, textvariable=self.effect_var, state="disabled",
                                      width=19, style="Filter.TCombobox")
        self.cb_effect['values'] = ["Все"]
        self.cb_effect.bind("<<ComboboxSelected>>", lambda e: self.apply_filters())
        self.cb_effect.pack(side=tk.LEFT)

        # 3. Таблица результатов
        table_container = tk.Frame(self, bg=self.COLOR_CARD, padx=2, pady=2,
                                   highlightthickness=1, highlightbackground=self.COLOR_BORDER)
        table_container.pack(fill=tk.BOTH, expand=True, padx=15, pady=6)

        columns = ("strain", "pos", "type", "ref", "alt", "effect", "stop", "aa_count")
        self.tree = ttk.Treeview(table_container, columns=columns, show="headings", selectmode="browse")

        self.tree.tag_configure("even_strain", background=self.COLOR_ROW_EVEN)
        self.tree.tag_configure("odd_strain", background=self.COLOR_ROW_ODD)

        headers = {
            "strain": "Штамм", "pos": "Позиция", "type": "Тип",
            "ref": "Реф", "alt": "Прим", "effect": "Эффект",
            "stop": "Новый стоп", "aa_count": "Затронуто АК"
        }
        col_widths = {
            "strain": 200, "pos": 75, "type": 95, "ref": 65,
            "alt": 65, "effect": 180, "stop": 105, "aa_count": 95
        }

        for col, title in headers.items():
            self.tree.heading(col, text=title)
            self.tree.column(col, width=col_widths[col], anchor=tk.CENTER if col != "strain" else tk.W)

        # Подключаем стилизованный скроллбар
        scroll_y = ttk.Scrollbar(table_container, orient=tk.VERTICAL,
                                 command=self.tree.yview, style="Vertical.TScrollbar")
        self.tree.configure(yscrollcommand=scroll_y.set)

        self.tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scroll_y.pack(side=tk.RIGHT, fill=tk.Y)

        # 4. Статус-бар
        status_frame = tk.Frame(self, bg=self.COLOR_CARD, padx=15, pady=6,
                                highlightthickness=1, highlightbackground=self.COLOR_BORDER)
        status_frame.pack(fill=tk.X, side=tk.BOTTOM, padx=15, pady=(0, 10))

        self.status_var = tk.StringVar(value="Готов к работе. Выберите CSV-файл со штаммами.")
        lbl_status = tk.Label(status_frame, textvariable=self.status_var, bg=self.COLOR_CARD,
                              fg=self.COLOR_MUTED, font=("Segoe UI", 9))
        lbl_status.pack(side=tk.LEFT)

    def select_file(self):
        path = filedialog.askopenfilename(
            title="Выберите CSV со штаммами",
            filetypes=[("CSV файлы", "*.csv"), ("Все файлы", "*.*")]
        )
        if path:
            self.input_file = path
            name = path.split("/")[-1].split("\\")[-1]
            self.lbl_file.config(text=name, fg=self.COLOR_TEXT, font=("Segoe UI", 9, "bold"))
            self.btn_run.config(state=tk.NORMAL)

    def run_analysis(self):
        self.all_results.clear()

        try:
            strains = load_from_csv(self.input_file)
            if not strains:
                messagebox.showwarning("Внимание", "Файл пуст или повреждён.")
                return

            types_set = set()
            effects_set = set()

            for strain_idx, item in enumerate(strains):
                name = item['strain_name']
                ref = item['ref_dna']
                sample = item['sample_dna']

                muts, fs_info = analyze_mutations(ref, sample, item['start_pos'])
                fs_map = {pos: (stop_pos, count) for pos, stop_pos, count in fs_info}

                for m in muts:
                    pos = m['Позиция']
                    stop_pos, affected_aa = fs_map.get(pos, ('-', '-'))
                    stop_display = str(stop_pos) if stop_pos is not None else 'Не найден'

                    types_set.add(m['Тип'])
                    effects_set.add(m['Эффект'])

                    self.all_results.append({
                        'strain_idx': strain_idx,
                        'Штамм': name,
                        'Эталон ДНК': ref,
                        'ДНК образца': sample,
                        'Позиция': pos,
                        'Тип': m['Тип'],
                        'Реф. аллель': m['Реф. аллель'],
                        'Прим. аллель': m['Прим. аллель'],
                        'Эффект': m['Эффект'],
                        'Новый стоп-кодон': stop_display,
                        'Затронуто аминокислот': affected_aa
                    })

            # Заполняем фильтры найденными значениями
            self.cb_type['values'] = ["Все"] + sorted(list(types_set))
            self.cb_type.current(0)
            self.cb_type.config(state="readonly")

            self.cb_effect['values'] = ["Все"] + sorted(list(effects_set))
            self.cb_effect.current(0)
            self.cb_effect.config(state="readonly")

            self.apply_filters()
            self.btn_export.config(state=tk.NORMAL)

        except Exception as err:
            messagebox.showerror("Ошибка обработки", str(err))

    def apply_filters(self):
        self.tree.delete(*self.tree.get_children())
        self.displayed_results.clear()

        selected_type = self.type_var.get()
        selected_effect = self.effect_var.get()

        for row in self.all_results:
            if selected_type != "Все" and row['Тип'] != selected_type:
                continue
            if selected_effect != "Все" and row['Эффект'] != selected_effect:
                continue

            self.displayed_results.append(row)
            row_tag = "even_strain" if row['strain_idx'] % 2 == 0 else "odd_strain"

            self.tree.insert("", tk.END, values=(
                row['Штамм'], row['Позиция'], row['Тип'], row['Реф. аллель'],
                row['Прим. аллель'], row['Эффект'], row['Новый стоп-кодон'],
                row['Затронуто аминокислот']
            ), tags=(row_tag,))

        # Пересчитываем статусную строку
        types = Counter(r['Тип'] for r in self.displayed_results)
        fs_count = sum(1 for r in self.displayed_results if r['Эффект'] == 'Сдвиг рамки считывания')
        self.status_var.set(
            f"Отображено: {len(self.displayed_results)} из {len(self.all_results)}   |   "
            f"Замен: {types.get('Замена', 0)}   |   "
            f"Делеций: {types.get('Делеция', 0)}   |   "
            f"Вставок: {types.get('Вставка', 0)}   |   "
            f"Сдвигов рамки: {fs_count}"
        )

    def export_csv(self):
        if not self.displayed_results:
            messagebox.showwarning("Экспорт", "Нет данных для сохранения.")
            return

        export_path = filedialog.asksaveasfilename(
            defaultextension=".csv",
            filetypes=[("CSV файл", "*.csv")],
            initialfile="mutations_report.csv",
            title="Сохранить отчет как"
        )
        if not export_path:
            return

        try:
            clean_results = [
                {k: v for k, v in row.items() if k != 'strain_idx'}
                for row in self.displayed_results
            ]

            # Экспортирует ровно то, что сейчас отфильтровано
            save_mutations_to_csv(clean_results, export_path)
            messagebox.showinfo("Готово",
                                f"Отчет успешно сохранен ({len(clean_results)} строк):\n{export_path}")
        except Exception as e:
            messagebox.showerror("Ошибка сохранения", str(e))


if __name__ == "__main__":
    app = DNAMutationApp()
    app.mainloop()