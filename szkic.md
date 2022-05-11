---
Title: Ideas for scripts to work with CCR
Tags: doktorat, programowanie, skrypt, python
Author: PinaMarzec
Source:
----


General idea for scripts to work with CCR experiments
==

Przygotowanie do pracy z CCR
-- 
1. Instrukcja
2. lista programów i bibliotek
     - python 2 i 3 
     - qMDD
     - sparky
     - 
 3. Pliki wejściowe
     - przesunięcia chemiczne badanego białka w formacie BMRB/NMRSTAR....
     - lista eksperymentów - experiment_set.txt
 4. 



# I. Puszczanie wielu eksperymentów na Brukerze (Python 2 / Jython)




# II. Przygotowywanie widm (Python 3 / Bash)
- pamiętać o komunikatach błędu!
- ścieżki do wyboru - czy FT czy CS, czy kopiowanie czy nie, itp
  
1. plik z informacjami gdzie są dane, gdzie je przenieść i które jest do którego CCR'a i jaka wersja:
     pre_experiments_set.txt
2. kopiowanie plików 
3. przygotowywanie plików qMDD (i ich edycja) 
     - Czy folder musi się nazywać coś.proc ???????  (bash wywala błąd przy robieniu czegoś z tym folderem - działa, ale nie wiem czy mogą być potem jakieś problemy - nie musi!)
     - proc.sh może mieć dwa zestawy danych, razem pracować nad wersją auto i cross
     - dwie wersje: 1 - FT, 2 - CS
     - tworzenie pliku experiments_set.txt (dane z pre_experiments_set, a resztę wpisuje program):
         - typ CCR
         - wersja CCR
         - ilość skanów dla obu wersji
         - relatywna pozycja pików/jąder dla jakiego aa, jaki kąt
4. puszczanie skryptu proc.sh - 
5. porządkowanie zbieranie widm razem



# III. Przygotowywanie list pików 
    
## IIIa. przygotowywanie startowych list pików (Python 2/3?) i plików .save i .proj
- Wszystko.py przygotowujące list pików w formacie Sparkiego, bez rozchwiania
- w oparciu o listę przesunięć chemicznych
- sprawdza jakie listy pików są potrzebne
- pytać czy są listy pików czy tworzymy nowe

## IIIb. pikowanie wielu widm (Python 3) - razem z powyższym
1. przygotowywanie plików .save i .proj 
   - pozycje pików z listy pików (chyba, że wolimy wczytywać przez Sparkiego)
   - synchronizacja wszystkich widm
   - nadanie odpowiednich nazw 
   - ustawianie ct

2. Skrypt do Sparkygo? (Python 2/3?)  - to może być trudne
    - otwieranie kilku widm na raz (fm)
    - wstawianie listy pików (rp) na widma 'a' - chyba, że pliki .save będą dobrze przygotowane
    - centrowanie pików (pc) w 'a'
    - zapisywanie listy pików z 'a'
    - przekopiowanie pików (pa oc op) na widma 'x' 
    - cetrowanie silnych pików w 'x'
    - zapisywanie listy pików z 'x'
    - zapisywanie projektu (js)
> porzucilyśmy ten pomysł i teraz mamy centrowanie za pomocą skryptu, który liczy poziom szumów i okresla czy dany pik jest widoczny (read_ucsf)

### read_ucsf
1. Wczytywanie danych
   
~~~
    read_ucsf script reading ucsf file and peak list (in Sparky format), check peak intensity of peak and if it is not in the highest position - move it. 
    After this script prints peak lists in ppm value and points value.
    
    For run script type in command line:
        python3 read_ucsf [ucsf path] [peak list path]
      
    additionaly you can add:
        --np [num]\t- to change number of points for calculate noise level: N^(spectra dimentionality + 1); normally is N = 10
        --plevel [num]\t- if you know level when starting appear, add this with scientific numer notation e.g. 1e+7
        --noRemove\t- add this if you do not want remove invisible peaks
        --onlypoints\t- add this if you want only change ppm value to points value
~~~

2. Wyznaczanie poziomu szumu dla każdego piku:
   - bierzemy przesunięcia chemiczne po protonie dla każdego piku
   - 
> nie, losujemy N^(spectra dimentionality+1)^ (normally is N = 10) i sprawdzamy czy ten punkt jest w pobliżu piku

   - zczytujemy wartości punktów z obszaru bez pików o danym przesunięciu protonu - w jakiejś odległości od innych pików
   - ~~automatyczny pikpiking~~
   - sprawdzić jak wyciągać dane o wartości punktów bezpośrenio z ucsf (przeczytać manuala Sparkiego)

3. Peak centering 
   - sprawdza wartości (intensywności) piku i przestrzeni o jeden w każdą strone
   - jeśli najwyzsze miejsce jest w innym miejscu niż pik to przesuwa pik i sprawdza ponownie (sprawdza maxymalnie o 2 punkty od orginalnego piku)
   - na koniec sprawdza czy to najwyzsza pozycja jesli nie to wraca do orgilanej pozycji
4. dodatkowe opcje:
  - [x] - wczytanie podanego przez użytkownika poziomu odcięcia
  - [ ] - ReportBox - plik ze wszystkimi infomacjami o przeprowadzonych obliczeniach
  - [ ] - ustawienia outputu przez użytkownika - nazwy i folderu
  - [ ] - wyswietlanie uzywanych dodatkowych opcji 


### read_point_list_compere - skrypt do porównywania pozycji pików na podstawie list pików 
   - [x] - porównywanie na podstawie pozycji z punktach
   - [x] - wczytywanie ścieżki dostepu do plików z pliku źródłowego 
   - [x] - wczytywanie i porównywanie kilku zestawów danych 
   - [ ] - wyświetlanie wysokości jeśli to możliwe
   - [x] - jeśli punkty są w postaci float to sprawdza czy różnica jest mniejsza niż pól punkta jeśli tak to "OK", jeśli nie to "Change", również w dodatkowej funkcji compere2list

## IIIc. przenoszenie danych między komputerami (Bash)
- konkretny schemat ułożenia danych:

   |- dowolna nazwa
       |-OrginFiles                    (folder z plikami z TopSpina)
           |- 1_CCR_1_a                nr folderu w TopSpinie_CCR_nr ccr>_wersja ccr
               |- 1_CCR_1_a.proc       (folder z plikami z MDD)
           |- 2_CCR_1_x
               |- 2_CCR_1_x.proc
           |- 3_CCR_2_a 
               |- 3_CCR_2_a.proc
           |- 4_CCR_2_x
               |- 4_CCR_2_x
       |-Spectra                      (folder z widmami)
           |-save
           |-peaklists

WERSJA 2

   |- dowolna nazwa

       |-OrginFiles                    (folder z plikami z TopSpina)
           |- CCR_1                    (folder z plikami z MDD i folderami z plikami TopSpina w wersi auto i cross)   
               |- 1_a                  (folder z danymi z Topspina dla wersji auto; nazwa: nr folderu w TopSpinie>_a
               |- 2_x                  (folder z danymi z Topspina dla wersji cross; nazwa: nr folderu w TopSpinie>_x
               |- FT_MDD

           |- CCR_1                      
               |- 1_a
               |- 2_x     
           |- 2_CCR_1_x
               |- 2_CCR_1_x.proc
           |- 3_CCR_2_a 
               |- 3_CCR_2_a.proc
           |- 4_CCR_2_x
               |- 4_CCR_2_x
       |-Spectra                      (folder z widmami)
           |-save
           |-peaklist


- plik z opisem jak używać skryptu

- (do przemyślenia) skrypt Sparkiego do zapisywania danych w naszym schemacie ułożenia danych???



IV. Liczenie stałych CCR (Python 3)
--

### IVa. Potrzebne pliki
  1. listy pików 
      - format listy pików Sparkiego
  2. experiments_set.txt 
      - kolejność jąder zgodna z listą pików
  3. sekwencja
      - FASTA
### IVb. Kod
- [ ] wczytywanie listy pików - nie tylko intensywności, ale też pozycje pików by określać przykrywanie.
- [ ] sprawdzanie przykrywania się pików ()
- [ ] przeliczanie wartości stałych CCR 
- [ ] eksport wyników 
- [ ] rysowanie wykresów
- [ ] jeśli są to dane dla bialka zwiniętego możemy porównac do teorii - kompatypilność z danymi teoretycznymi i ewentulanie ze skryptem Ani do liczenia kątów
1. Klasy:
    - [x] CSpectrum   (informacje o widmach)
    - [x] CPeak       (informacje o pikach)
        - [ ] PrintPeak - metoda/funkcja do wypisywania informacji o konkretnym piku
    - [x] CSequence   (informacje o sekwencji i wartości stałych CCR)
2. Funkcje:
    - [x] one2three - zmiana nazw aminikowasów z 1-literowego a 3-literowy
    - [x] three2one - zmiana nazw aminikowasów z 3-literowego a 2-literowy
    - 
    - [ ] ReadInputFile - wczytywanie plików wejściowych do klas CSpectrum, CPeak, CSequence 
        - [x] ReadExpSet - wczytywanie danych z pliku experiment_set.txt do klasy CSpectrum
        - [x] ReadSequnce - wczytywanie sekwencji aminokwaów do klasy CSequence --> na razie tylko z FASTA
        - [x] ReadPeakList - wczytywanie pozycji pików oraz ich wysokości do klas CPeak i CSpectrum
    - [ ] CheckOverlap - porównywanie odległości między pikamy w danych wymiarach i uzupełaniaie informacji w klasie CPeak
        - jeśli odległość w danych wymiarach będzie mniejsza niż .... to oznaczamy jako przykrywanie
        - oznaczamy od razu w obu piki 
    - [ ] CalcCCRRate - obliczanie stałych CCR 
    - [ ] WriteCCRRate - wypisywanie pliku ze stałymi CCR i adnotacjiami jakościowymi (przykrywanie, intensywność pików)
    - [ ] Write... - wypisywanie pośrednich plików, np. informacje o przykrywaniach pików w danych widmach
IVc. Pliki wyjściowe
    - CCRRate.txt lub .csv - sekwencja, stałe CCR, podstawowe adnotacje (w rozszerzeniu csv lub txt)
    - log.txt - plik ze wszsytkim co wypisuje program w terminalu
    - .... - plik z pikami dla każdego widma i zaznaczonymi wszystkimi przykrywaniami









