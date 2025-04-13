# Установка необходимых пакетов
install.packages("citation")  # Установка пакета citation
install.packages("tidyverse") # Установка tidyverse для работы с данными
install.packages("sf")        # Установка sf для работы с пространственными данными
install.packages("mapview")   # Установка mapview для визуализации карт
install.packages("tmap")      # Установка tmap для создания тематических карт
install.packages("terra")     # Установка terra для работы с растровыми данными

# Подгрузка необходимых библиотек
library(tidyverse)    # Для работы с данными
library(sf)           # Для работы с пространственными данными
library(mapview)      # Для визуализации карт
mapviewOptions(fgb = FALSE)  # Отключение формата FGB
mapviewOptions(viewer.suppress = FALSE)  # Разрешение отображения карты
library(tmap)         # Для создания тематических карт
library(rcartocolor)  # Для использования цветовых палитр
library(rmapshaper)   # Для обработки пространственных данных
library(terra)        # Для работы с растровыми данными
library(styler)       # Для форматирования кода
 

# ===================
# 1. Подготовка данных
# Загрузка растра населения
rast <- rast("Data/Raw/GHS_POP_E2020_GLOBE_R2023A_54009_100_V1_0_R4_C22.tif") |>
  aggregate(fact = 10, fun = "sum", na.rm = TRUE)  # Агрегация растра

# Проецирование растра в равновеликую проекцию Альберса
albers <- crs("+proj=aea +lat_1=44 +lat_2=48 +lat_0=46 +lon_0=39 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
raster_albers <- project(rast, y = albers, res = 1000, method = "average")  # Проекция растра

# Проверка суммарного населения
global(rast, "sum", na.rm = TRUE)

# Загрузка границ Краснодарского края и Республики Адыгея
krasnodar <- st_read("Data/Raw/Субъекты_России.shp") |>
  filter(region %in% c("Краснодарский край", "Республика Адыгея")) |>
  st_transform(crs = st_crs(albers))  # Преобразование в проекцию Альберса

# Обрезка растра по границам Краснодарского края и маскирование
rast_mask <- raster_albers |>
  terra::crop(vect(krasnodar)) |>
  mask(vect(krasnodar))  # Обрезка и маскировка растра

# Проверка суммарного населения
global(rast_mask, "sum", na.rm = TRUE)


mapview(rast_mask)  # Визуализация растра

# ===================
# 2. Классификация данных
# Матрица классификации (урбанизированные районы)
rcl <- matrix(
  c(-Inf, 1500, 0,  # Население < 1500 -> 0
    1500, Inf, 1),  # Население >= 1500 -> 1
  ncol = 3, byrow = TRUE
)

# Классификация растра
classified_rast <- classify(rast_mask, rcl)

# Визуализация результатов классификации
plot(classified_rast)

# Преобразование классифицированного растра в полигоны
rast_polys <- as.polygons(classified_rast, na.rm = TRUE, aggregate = TRUE)

# Преобразование SpatVector в объект sf
rast_sf <- st_as_sf(rast_polys) |>
  st_cast("POLYGON") |>
  filter(GHS_POP_E2020_GLOBE_R2023A_54009_100_V1_0_R4_C22 == 1)  # Фильтрация полигонов

mapview(rast_sf)  # Визуализация полигонов

# Расчет населения для каждого полигона
union <- terra::extract(rast_mask, vect(rast_sf), fun = sum, na.rm = TRUE)

# Добавление информации о населении к слою с полигонами
total_pop <- rast_sf |>
  mutate(total_pop = union[, 2]) |>
  filter(total_pop > 50000)  # Фильтрация полигонов с населением > 50 тыс.

# ===================
# 3. Заполнение пропусков в полигоне
# Функция для заполнения пропусков
fill_gaps <- function(test_polygon, rast_mask) {
  # Растеризация полигона
  test_raster <- rasterize(vect(test_polygon), rast_mask)
  test_raster[is.na(test_raster)] <- 0
  
  # Расчет суммы в окне 3x3
  test_focal <- focal(test_raster, w = matrix(1, 3, 3), fun = sum)
  
  # Классификация
  rcl <- matrix(
    c(-Inf, 4, 0,  # Сумма <= 4 -> 0
      4, Inf, 1),  # Сумма > 4 -> 1
    ncol = 3, byrow = TRUE
  )
  classified_focal <- classify(test_focal, rcl)
  
  # Объединение растров
  res_raster <- test_raster + classified_focal
  res_raster[res_raster > 0] <- 1
  
  # Преобразование в полигон
  res_polygons <- res_raster |>
    as.polygons(na.rm = TRUE, aggregate = TRUE) |>
    st_as_sf() |>
    filter(layer == 1)
  
  return(res_polygons)
}

# Итерационное заполнение пропусков
itog_sloy <- tibble()  # Итоговый слой
for (f in 1:nrow(total_pop)) {
  old_poly <- total_pop[f, ]
  new_poly <- fill_gaps(old_poly, rast_mask)
  
  i <- 1
  while (!identical(st_area(new_poly), st_area(old_poly))) {
    print(i)
    old_poly <- new_poly
    new_poly <- fill_gaps(old_poly, rast_mask)
    i <- i + 1
  }
  
  itog_sloy <- rbind(itog_sloy, new_poly)
}

# Удаление пересечений между полигонами
urban_centers <- st_difference(itog_sloy)

# ===================
# 4. Классификация полугородских районов
# Матрица классификации (полугородские районы)
rcl_suburbia <- matrix(
  c(-Inf, 300, 0,  # Население < 300 -> 0
    300, Inf, 1),  # Население >= 300 -> 1
  ncol = 3, byrow = TRUE
)

# Классификация растра
classified_rast_suburbia <- classify(rast_mask, rcl_suburbia)
mapview(classified_rast_suburbia)

# Преобразование в полигоны
rast_polys_suburbia <- as.polygons(patches(classified_rast_suburbia, directions = 8, zeroAsNA = TRUE), na.rm = TRUE, aggregate = TRUE)

rast_sf_suburbia <- st_as_sf(rast_polys_suburbia) |>
  st_cast("MULTIPOLYGON")

# Расчет населения для каждого полигона
union_suburbia <- terra::extract(rast_mask, vect(rast_sf_suburbia), fun = sum, na.rm = TRUE)

rast_sf_suburbia |> mutate(total_pop = union_suburbia[, 2]) -> total_pop_suburbia
# Добавление информации о населении к слою с полигонами
total_pop_suburbia <- filter(total_pop_suburbia, total_pop >= 5000)  # Фильтрация полигонов с населением > 5 тыс.
mapview(total_pop_suburbia)

# Вычитание урбанизированных районов
suburbia_sloy <- rmapshaper::ms_erase(total_pop_suburbia, urban_centers) |>
  st_cast("MULTIPOLYGON") |>
  st_cast("POLYGON")
mapview(suburbia_sloy)
# ===================
# 5. Классификация сельской местности
# Матрица классификации (сельская местность)
rcl_rural <- matrix(
  c(-Inf, 0, 0,  # Население == 0 -> 0
    0, Inf, 0),  # Население > 0 -> 1
  ncol = 3, byrow = TRUE
)
# Классификация растра
classified_rast_rural <- classify(rast_mask, rcl_rural)


# Преобразование в полигоны
rast_polys_rural <- as.polygons(classified_rast_rural, na.rm = TRUE, aggregate = TRUE)
rast_sf_rural <- st_as_sf(rast_polys_rural) |>
  st_cast("POLYGON")

# Вычитание урбанизированных районов и полугородских районов
rural_sloy <- rmapshaper::ms_erase(rast_sf_rural, urban_centers) |>
  rmapshaper::ms_erase(suburbia_sloy) |>
  st_cast("MULTIPOLYGON") |>
  st_cast("POLYGON")
mapview(rural_sloy) +
  mapview(urban_centers)



# ===================
# 6. Объединение слоев
OCED_classification <- rbind(
  urban_centers |> dplyr::select(geometry) |> mutate(type = "urban center"),
  suburbia_sloy |> dplyr::select(geometry) |> mutate(type = "urban cluster"),
  rural_sloy |> dplyr::select(geometry) |> mutate(type = "rural")
) |>
  mutate(type = factor(type, levels = c("urban center", "urban cluster", "rural")))

# Визуализация итоговой классификации
tm_shape(OCED_classification) +
  tm_polygons("type")
mapview(OCED_classification)
# Сохранение итогового слоя
write_sf(OCED_classification, "Data/Ready/OECD_classification.shp")


# ===================
# 7. Расчет численности населения по типам
# Создание нового слоя с пересечением полигонов типов и полигонов регионов
regions_types <- st_intersection(krasnodar, OCED_classification) |>
  dplyr::select(region, fo, type)

# Посчитали по населению полигоны через растр
regions_types_population <- terra::extract(rast_mask, vect(regions_types), fun = sum, na.rm = TRUE)

# Добавили колонку с населением в слой
regions_types <- regions_types |>
  mutate(population = regions_types_population[, 2])

# Группировка и агрегирование данных
types_population_for_regions <- regions_types |>
  st_drop_geometry() |>
  group_by(region, type) |>
  summarise(population = sum(population, na.rm = TRUE)) |>
  ungroup()

# Добавим колонку с долей населения каждого типа
types_population_for_regions <- types_population_for_regions |>
  group_by(region) |>
  mutate(proportion_population = round(population / sum(population) * 100, 2)) |>
  ungroup()

# Сделаем график
ggplot(types_population_for_regions, aes(x = region, y = proportion_population, fill = type)) +
  geom_col(position = "stack") +
  labs(title = "Доля населения по типам региона",
       x = "Регион",
       y = "Доля населения (%)",
       fill = "Тип региона")


