# -*- coding: utf-8 -*-

"""
Created on Mon Nov  8 10:39:38 2021

@author: VNG
"""

import os
import glob
import math as m
import configparser
import random
import platform
import subprocess
from datetime import datetime
from string import ascii_letters
from zipfile import is_zipfile
from shutil import rmtree

import pandas as pd
from sklearn.metrics.pairwise import haversine_distances


# ==================================================================
# Общие константы, которые не предполагается менять пользователю
# ==================================================================
# радиус Земли
EARTH_KM = 6371
# сколько км в 1 градусе меридианов (зависит от того, каким мы считаем радиус Земли)
ONE_DEGREE_LONGITUDE = EARTH_KM / 180 * m.pi  # 111.1949 км при радиусе = 6371
COEFF = 1.01


# =================== Основной класс ===============================
class FindTopClosest:
    """
    Класс для нахождения ближайших точек по геокоординатам.
    На вход:
        архив в ZIP с текстовым файлом настроек и 2 Excel с геокоординатами
        1. точек интереса; 2. точек, до которых хотим дотянуться.
    Результат:
        также в ZIP-архиве: 1.Excel-список ближайших нйденных с учетом получ.
        на входе параметров (расстояние, кол-во макс.найд.); 2.текстовый файл
        со значениями макс.возможного кол-ва найденных для каждой т.интереса
        при заданном расстоянии.
    """

    def __init__(self,
                 in_file,
                 target_folder):
        """
        Инициализация представителя класса"""
        self.in_file = os.path.abspath(in_file)
        self.start_time = datetime.now()
        print(f'Время начала работы: {self.start_time}')
        datetime_txt = self.start_time.strftime('%Y%m%d_%H%M%S')
        folder_out = os.path.join(target_folder,
                                  "geo_dist_" + datetime_txt)
        self.folder_out = folder_out    # папка для сохранения результатов
        # создаем папку для сохранения результатов.
        create_folder(folder_out)

        # Временные файлы
        self.temp_path = os.path.join(folder_out, 'temp')
        create_folder(self.temp_path)

        # пути к рабочим файлам (читаем из конф.файла след.шагом)
        self.fname_POI = None
        self.fname_check = None
        # параметры для поиска (также из конф.файла чуть позже)
        # радиус поиска объектов
        self.interest_radius_km = None
        # кол-во ближайших по расстоянию, которых возвращаем как результат
        self.top_n = None

        # Результат работы в заархивированном виде -- полное имя файла
        self.arch_result = None
        self.status = 'init'

    def read_config(self):
        """
        Распаковываем архив и читаем настройки конфигурации из
        конфигурационного файла.
        """
        # проверяем, а zip-архив ли нам передали
        if not is_zipfile(self.in_file):
            print(f'Не является архивом ZIP или неверный путь к файлу: {self.in_file}.')
            return
        # разархивируем во временную папку
        unzip_exit_stat = unzip_unzip(self.in_file, self.temp_path)
        if unzip_exit_stat.returncode > 0:
            return

        config = configparser.ConfigParser()
        config.read(os.path.join(self.temp_path,
                                 'geo_dist_config.ini'), encoding='utf-8')

        # Читаем из файла конфигурации имена Excel-файлов с данными...
        self.fname_POI = os.path.join(self.temp_path,
                                      config["Paths_in"]["f_POI_name"])
        self.fname_check = os.path.join(self.temp_path,
                                        config["Paths_in"]["f_check_pts_name"])
        # ... и параметры обработки
        # радиус интереса в км
        self.interest_radius_km = float(config["Params"]["interest_radius_km"])
        # кол-во ближайших по расстоянию объектов
        self.top_n = int(config["Params"]["top_n"])
        self.status = 'read_config'

    def calculate(self):
        """
        Проводим основные вычисления, формируем итоговые файлы.
        """
        # Вычисляем дельту градусов по широте (меридианы считаем одинаковыми)
        #
        EPS_LATTITUDE = COEFF * self.interest_radius_km / ONE_DEGREE_LONGITUDE
        # Пути к файлам
        df_poi = pd.read_excel(self.fname_POI, engine='openpyxl')
        df_check = pd.read_excel(self.fname_check, engine='openpyxl')
        # Исправляем долготу, если вдруг есть Западная (с минусом)
        df_poi = fix_long360(df_poi)
        df_check = fix_long360(df_check)

        # Удаляем временные файлы
        clear_tmp(self.temp_path)
        # Удаляем временную папку
        rmtree(self.temp_path, ignore_errors=True)

        df_fin = pd.DataFrame()

        out_max_qty_file = os.path.join(self.folder_out,
                                   txt_tstamp(f'r={self.interest_radius_km}_max_quantity.txt'))
        base_fname_check = os.path.basename(self.fname_check)
        write_txt_a(out_max_qty_file, f'Кол-во точек из файла {base_fname_check} попало в рассмотрение:\n', mode='w')

        for row in df_poi.iloc[:, :3].itertuples():
            # Берем очередную строку с данными о точке интереса (POI = Point Of Interest)
            __, poi_name, poi_latt, poi_long = row
            poi_radians = [m.radians(poi_latt), m.radians(poi_long)]

            # == Фильтруем для нее выборку с точками заданного квадрата
            # Вычисляем дельту по долготе (дельта по долготе вычисл. с помощью
            # длины параллели, а длина окружности параллелей уменьш. к полюсам)
            eps_longitude =  COEFF * self.interest_radius_km * one_km_to_degree(poi_latt)
            # После применения фильтра выводим размеры получившегося датафрейма
            vsp_local = filter_df(df_check, poi_latt, poi_long, EPS_LATTITUDE, eps_longitude)
            t_out = f'Для {poi_name} было отфильтровано {vsp_local.shape[0]}.\n'
            write_txt_a(out_max_qty_file, t_out)

            # Вычисляем расстояния в словарь, проходясь по всем отфильтрованным значениям координат
            dist_calculated = {}
            # Если ничего не попало с заданными условиями, то
            # записываем одну строчку - "Не найдено" и р-ние 99999 км.
            if vsp_local.shape[0] == 0:
                df_calc = pd.DataFrame({'id_check_point':['not_found'],
                                        'distance_km': [99999]})
            else:
                for row_vsp in vsp_local.itertuples():
                    __, vsp_name, vsp_latt, vsp_long = row_vsp
                    vsp_radians = [m.radians(vsp_latt), m.radians(vsp_long)]
                    distance = EARTH_KM * haversine_distances([poi_radians, vsp_radians])[0, 1]
                    dist_calculated[vsp_name] = distance

                # Составляем датафрейм
                df_calc = pd.DataFrame(data=dist_calculated.values(), index=dist_calculated.keys())
                df_calc.sort_values(by=[0], ascending=True, inplace=True)
                df_calc.reset_index(inplace=True)
                df_calc.rename(columns={'index': 'id_check_point', 0: 'distance_km'}, inplace=True)

                # Убираем тех, кто попал с превышением расстояния
                # (из углов предварительного квадрата, например)
                df_calc = df_calc[df_calc['distance_km'] <= self.interest_radius_km]
                # Если после отбора по расстоянию отбросили все найденные до этого, то создаем строку
                if df_calc.shape[0] == 0:
                    df_calc = pd.DataFrame({'id_check_point':['not_found'],
                                        'distance_km': [99999]})
                # Выделяем лидирующие позиции
                df_calc = df_calc.head(self.top_n)
            df_calc['poi_name'] = poi_name
            # Дописываем в конец датафрейм с результатами
            df_fin = df_fin.append(df_calc)

        res_xls = os.path.join(self.folder_out,
                               txt_tstamp(f'Расчет_{self.interest_radius_km}км_{self.top_n}макс.xlsx'))
        df_fin.to_excel(res_xls)
        self.status = 'calculate'
        print(f'Работа завершена: {datetime.now()}.\nПродолжительность работы: {datetime.now()-self.start_time}')

# ============== MAIN для сервера, по образцу предыдущих проектов ===========
class geo_dist_serv:
    """
    Класс для работы модуля на сервере.
    """
    def __init__(self, file_name, dest_folder):
        self.file_name = file_name
        self.dest_folder = dest_folder

    def start(self):
        """
        Основная функция, собирающая вместе все шаги по обработке,
        вычислениям и записи в файлы.
        """

        find_closest = FindTopClosest(self.file_name, self.dest_folder)
        find_closest.read_config()
        find_closest.calculate()


# ========== вспомогательные функции ========================================
def create_folder(new_folder):
    """
    Создаем папку
    ver ocr
    """
    try: # пытаемся создать всю структуру папок сразу
        os.makedirs(new_folder)
    except (OSError, PermissionError):
        pass

def generate_id(str_len=5):
    """
    Генерит случайную последовательность из [a-zA-Z]
    длиной str_len.
    ver ocr
    """
    random.seed()
    return "".join(random.sample(ascii_letters, str_len))

def unzip_unzip(zip_filename, extract_folder):
    """
    (вариант 2)
    С помощью unzip
    извлекает файлы из архива zip_filename
    в указанную папку extract_folder
    ver ocr
    """
    params = ['unzip', '-qq', zip_filename, '-d', extract_folder]
    return subprocess.run(params, check=True)

def clear_tmp(tmp_path):
    """
    Очистка папки с файлами
    """
    for path in glob.glob(os.path.join(tmp_path, '*')):
        try:
            os.remove(path)
        except:
            pass

def out_file(f_name, suffix='_geo_dist'):
    """Формируем имя для файла результата"""
    before = os.path.splitext(f_name)
    after = before[0] + suffix + before[1]
    return after

def pack_files(arch_name, source_folder):
    """
    Переходим в папку, содержащую файлы для архивации,
    архивируем с помощью zip и
    возвращаемся в ту папку, откуда переходили изначально.
    ver ocr
    """
    # формируем список из которого склеим через пробел одну длинную команду
    params = ['pushd', source_folder, '&&',
              'zip', '-r', arch_name, '.', '&&',
              'popd']
    command = " ".join(params)
    print(command)
    # в зависимости от ОС разные терминалы. указывать явно нужно только для Linux,
    # т.к. там по умолчанию используется /bin/sh, в котором нет pushd
    if platform.system() != 'Windows':
        ret = subprocess.run([command], shell=True, executable='/bin/bash', check=True)
    else:
        ret = subprocess.run(params, shell=True, check=True)
    return ret


def one_km_to_degree(lattitiude_deg):
    """
    Возвращает, сколько градусов в 1 км для заданной широты.
    Если передали широту больше 90, то возвращаем None
    """
    if lattitiude_deg > 90:
        return None

    one_degree_latt = round(ONE_DEGREE_LONGITUDE * m.cos(lattitiude_deg * m.pi / 180), 8)
    if one_degree_latt > 0:
        return 1 / one_degree_latt

    return 0

def filter_df(df, lattit, longit, eps_latit, eps_longit):
    """
    Фильтруем датафрейм с координатами по заданному значению координат (longit, lattit)
    заданной погрешностью:
        eps_latit по широте;
        eps_longit по долготе.
    """
    df_filtered = df[(abs(df['Широта'] - lattit) <= eps_latit) &
                     (abs(df['Долгота'] - longit) <= eps_longit)
                    ]
    return df_filtered

def fix_long360(df, longit='Долгота'):
    """
    Для заданного датафрейма поправляем столбец с долготой.
    Если есть отриц.значения, то прибавляем к ним 360, чтобы можно было сравнивать с соседями.
    """
    df[longit] = df[longit].apply(lambda x: x+360 if x < 0 else x)
    return df

def write_txt_a(f_name: str, text, mode='a'):
    """Дописываем в текстовый файл список строк"""
    with open(f_name, mode, encoding='utf8') as f_out:
        f_out.writelines(text)

def txt_tstamp(filename: str, cur_time=None) -> str:
    """
    22.10.2021 vng
    В имя файла filename добавляет в конец (перед расширением) штамп времени
    в формате %Y-%m-%d_%H%M%S. По умолчанию -- текущее время,
    можно также передать заданное в виде datetime(2020, 12, 15, 7, 10, 45)
    """
    if cur_time is None:
        cur_time = datetime.now()
    curtime_txt = cur_time.strftime('%Y-%m-%d_%H%M%S')
    begining, ending = os.path.splitext(filename)
    new_filename = begining + "_" + curtime_txt + ending
    return new_filename
