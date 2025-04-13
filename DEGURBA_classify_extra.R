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
rast_polys <- as.polygons(patches(classified_rast, direction=4, zeroAsNA = TRUE), na.rm = TRUE, dissolve = TRUE)

# Преобразование SpatVector в объект sf
rast_sf <- st_as_sf(rast_polys)

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

mapview(urban_centers)







# Матрица классификации для плотных городских кластеров
rcl_dense_cluster <- matrix(
  c(-Inf, 1500, 0,  # Плотность < 1500 -> 0
    1500, Inf, 1),  # Плотность ≥ 1500 -> 1
  ncol = 3, byrow = TRUE
)

# Классификация растра
classified_rast_dense_cluster <- classify(rast_mask, rcl_dense_cluster)

# Преобразование в полигоны
dense_cluster_patches <- as.polygons(patches(classified_rast_dense_cluster, directions = 4, zeroAsNA = TRUE), na.rm = TRUE, dissolve = TRUE)
dense_cluster_sf <- st_as_sf(dense_cluster_patches)

# Расчет населения для каждого полигона
pop_dense_cluster <- terra::extract(rast_mask, vect(dense_cluster_sf), fun = sum, na.rm = TRUE)
dense_cluster_sf <- dense_cluster_sf |>
  mutate(total_pop = pop_dense_cluster[, 2]) |>
  filter(total_pop >= 5000 & total_pop < 50000)  # Фильтрация по населению [5000, 50000)


# Итерационное заполнение пропусков
dense_clusters_sloy <- tibble()  # Итоговый слой
for (f in 1:nrow(dense_cluster_sf)) {
  print(f)
  old_poly <- dense_cluster_sf[f, ]
  new_poly <- fill_gaps(old_poly, rast_mask)
  
  i <- 1
  while (!identical(st_area(new_poly), st_area(old_poly))) {
    print(i)
    old_poly <- new_poly
    new_poly <- fill_gaps(old_poly, rast_mask)
    i <- i + 1
  }
  
  dense_clusters_sloy <- rbind(dense_clusters_sloy, new_poly)
}


mapview(dense_cluster_sf)+
  mapview(dense_clusters_sloy)









# Матрица классификации для полуплотных городских кластеров
rcl_semi_dense_cluster <- matrix(
  c(-Inf, 900, 0,  # Плотность < 900 -> 0
    900, Inf, 1),  # Плотность ≥ 900 -> 1
  ncol = 3, byrow = TRUE
)

# Классификация растра
classified_rast_semi_dense_cluster <- classify(rast_mask, rcl_semi_dense_cluster)

# Преобразование в полигоны
semi_dense_cluster_patches <- as.polygons(patches(classified_rast_semi_dense_cluster, directions = 4, zeroAsNA = TRUE), na.rm = TRUE, dissolve = TRUE)
semi_dense_cluster_sf <- st_as_sf(semi_dense_cluster_patches)

# Расчет населения для каждого полигона
union_semi_dense_cluster <- terra::extract(rast_mask, vect(semi_dense_cluster_sf), fun = sum, na.rm = TRUE)
semi_dense_cluster_sf <- semi_dense_cluster_sf |>
  mutate(total_pop = union_semi_dense_cluster[, 2]) |>
  filter(total_pop >= 2500)  # Фильтрация по населению ≥ 2500

# Удаление пересечений с городскими центрами и плотными кластерами
semi_dense_cluster_sf <- semi_dense_cluster_sf |>
  rmapshaper::ms_erase(urban_centers) |>
  rmapshaper::ms_erase(dense_clusters_sloy)
mapview(semi_dense_cluster_sf)

# Удаление объектов в пределах 2 км от городских центров и плотных кластеров
buffer_2km <- st_buffer(rbind(urban_centers, dense_clusters_sloy), dist = 2000)
semi_dense_cluster_sf <- semi_dense_cluster_sf |>
  rmapshaper::ms_erase(buffer_2km)
mapview(semi_dense_cluster_sf)





# # Создание буфера в 2 км вокруг городских центров и плотных кластеров
# buffer_2km_urban <- st_buffer(st_union(rbind(urban_center_sf, dense_cluster_sf)), dist = 2000)
# 
# # Определение пригородных ячеек (пересечение оставшихся городских ячеек с буфером)
# suburban_sf <- st_intersection(remaining_urban, buffer_2km_urban)
# 
# # Удаление пустых геометрий
# suburban_sf <- suburban_sf[!st_is_empty(suburban_sf), ]
# 
# mapview(suburban_sf)


# Матрица классификации для сельских кластеров
rcl_rural_cluster <- matrix(
  c(-Inf, 300, 0,  # Плотность < 300 -> 0
    300, Inf, 1),  # Плотность ≥ 300 -> 1
  ncol = 3, byrow = TRUE
)

# Классификация растра
classified_rast_rural_cluster <- classify(rast_mask, rcl_rural_cluster)

# Преобразование в полигоны
rural_cluster_patches <- as.polygons(patches(classified_rast_rural_cluster, directions = 8, zeroAsNA = TRUE), na.rm = TRUE, dissolve = TRUE)
rural_cluster_sf <- st_as_sf(rural_cluster_patches)

# Создадим слой с пригородами
suburban <- st_filter(rural_cluster_sf, buffer_2km, .predicate = st_intersects)
mapview(suburban) +
  mapview(semi_dense_cluster_sf)

suburban <- suburban |>
  rmapshaper::ms_erase(urban_centers) |>
  rmapshaper::ms_erase(dense_clusters_sloy) |> 
  rmapshaper::ms_erase(semi_dense_cluster_sf)

mapview(semi_dense_cluster_sf)

# Расчет населения для каждого полигона
union_rural_cluster <- terra::extract(rast_mask, vect(rural_cluster_sf), fun = sum, na.rm = TRUE)
rural_cluster_sf <- rural_cluster_sf |>
  mutate(total_pop = union_rural_cluster[, 2]) |>
  filter(total_pop >= 500 & total_pop <= 4999) |>
  rmapshaper::ms_erase(suburban)











# Матрица классификации для сельских ячеек с низкой плотностью
rcl_low_density_rural <- matrix(
  c(-Inf, 50, 0,  # Плотность < 50 -> 0
    50, Inf, 1),  # Плотность ≥ 50 -> 1
  ncol = 3, byrow = TRUE
)

# Классификация растра
classified_rast_low_density_rural <- classify(rast_mask, rcl_low_density_rural)

# Преобразование в полигоны
low_density_rural_patches <- as.polygons(patches(classified_rast_low_density_rural, directions = 8, zeroAsNA = TRUE), na.rm = TRUE, dissolve = TRUE)
low_density_rural_patches_sf <- st_as_sf(low_density_rural_patches)


low_density_rural <- low_density_rural_patches_sf |>
  rmapshaper::ms_erase(urban_centers) |>
  rmapshaper::ms_erase(dense_clusters_sloy) |> 
  rmapshaper::ms_erase(semi_dense_cluster_sf) |> 
  rmapshaper::ms_erase(suburban) |> 
  rmapshaper::ms_erase(rural_cluster_sf)



  



# Матрица классификации для сельских районов с очень низкой плотностью
rcl_very_low_density_rural <- matrix(
  c(0, 50, 1,  # Население == 0 -> 0
    50, Inf, 0),  # Население > 0 -> 1
  ncol = 3, byrow = TRUE
)

# Классификация растра
classified_rast_very_low_density_rural <- classify(rast_mask, rcl_very_low_density_rural)


# Преобразование в полигоны
very_low_density_rural_patches <- as.polygons(patches(classified_rast_very_low_density_rural, directions = 8, zeroAsNA = TRUE), na.rm = TRUE, dissolve = TRUE)
very_low_density_rural_patches_sf <- st_as_sf(very_low_density_rural_patches)


very_low_density_rural <- very_low_density_rural_patches_sf |>
  rmapshaper::ms_erase(urban_centers) |>
  rmapshaper::ms_erase(dense_clusters_sloy) |> 
  rmapshaper::ms_erase(semi_dense_cluster_sf) |> 
  rmapshaper::ms_erase(suburban) |> 
  rmapshaper::ms_erase(rural_cluster_sf) |> 
  rmapshaper::ms_erase(low_density_rural)
  




# Проверка результата
mapview(semi_dense_cluster_sf, col.regions = "darkred")+
  mapview(rural_cluster_sf, col.regions = "darkgreen")+
  mapview(urban_centers, col.regions = "red")+
  mapview(suburban, col.regions = "yellow")+
  mapview(dense_clusters_sloy, col.regions = "gold4")+
  mapview(low_density_rural, col.regions = "palegreen3")+
  mapview(very_low_density_rural, col.regions = "palegreen")




# Стандартизация структуры данных для всех слоев
standardize_layer <- function(layer, type) {
  layer |>
    st_drop_geometry() |>
    select(total_pop = total_pop, everything()) |>
    mutate(type = type) |>
    st_as_sf(geometry = st_geometry(layer))
}

# Применение функции standardize_layer к каждому слою
urban_center_sf <- standardize_layer(urban_center_sf, "Urban Center")
dense_cluster_sf <- standardize_layer(dense_cluster_sf, "Dense Urban Cluster")
semi_dense_cluster_sf <- standardize_layer(semi_dense_cluster_sf, "Semi-Dense Urban Cluster")
suburban_sf <- standardize_layer(suburban_sf, "Suburban")
rural_cluster_sf <- standardize_layer(rural_cluster_sf, "Rural Cluster")
low_density_rural_sf <- standardize_layer(low_density_rural_sf, "Low Density Rural")
very_low_density_rural_sf <- standardize_layer(very_low_density_rural_sf, "Very Low Density Rural")

# Объединение всех типов
classification <- rbind(
  urban_center_sf,
  dense_cluster_sf,
  semi_dense_cluster_sf,
  suburban_sf,
  rural_cluster_sf,
  low_density_rural_sf,
  very_low_density_rural_sf
)

# Проверка результата
print(summary(classification))

# Визуализация
tm_shape(classification) +
  tm_polygons("type")