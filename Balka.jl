tcf(filename::String) = (@__DIR__) * "\\" * filename;
using Gmsh
using Gridap
using Gridap.Geometry
using Gridap.ReferenceFEs
using GridapGmsh

gmsh.initialize()

gmsh.model.add("model")

# Load a STEP file (using `importShapes' instead of `merge' allows to directly
# retrieve the tags of the highest dimensional imported entities):
path = "Balka.stp" |> tcf
v = gmsh.model.occ.importShapes(path)


# xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.occ.getBoundingBox(
#     v[1][1], v[2][2])

#     zmin=-13.0
 gmsh.model.occ.synchronize()
for plane in gmsh.model.getEntities(2)
    println(gmsh.model.getColor(plane...))
    gmsh.model.addPhysicalGroup(2, [plane[2]], -1, "$(plane[2])");

end

colors2 = [gmsh.model.getColor(plane...) for plane in gmsh.model.getEntities(2)] #1-15 и 3-13
d2 = gmsh.model.getEntities(2)

#Лево верх
bb_UpLeft = gmsh.model.getBoundingBox(d2[1]...) #верх лево
entitiesUp1_on_left_pad = gmsh.model.getEntitiesInBoundingBox(bb_UpLeft..., 1)
tagsUp1_on_left_pad=[et[2] for et in entitiesUp1_on_left_pad[1:5]] #Box для низа и верха совпадает, но почему то если брать бокс от верхней плоскости, то в неё входит 5 линий,
# а в нижний бокс всего 4. Если брать бокс от нижней плоскости то уже 4 линии от векрхней и 5 от нижней.
# точки лево верх
entities0_on_left_pad = gmsh.model.getEntitiesInBoundingBox(bb_UpLeft..., 0) # а с точками всё норм. Баг или фитча? 
tagsUp0_on_left_pad=[et[2] for et in entities0_on_left_pad[1:5]] #индексы точек левой верхней панели
#создаем физические группы с одинаковым названием для левой верхней грани
gmsh.model.addPhysicalGroup(2, [d2[1][2]], -1, "UpLeft")
gmsh.model.addPhysicalGroup(1, tagsUp1_on_left_pad, -1, "UpLeft")
gmsh.model.addPhysicalGroup(0, tagsUp0_on_left_pad, -1, "UpLeft")


#лево низ
bb_DownLeft = gmsh.model.getBoundingBox(d2[15]...)
entitiesDown1_on_left_pad = gmsh.model.getEntitiesInBoundingBox(bb_DownLeft..., 1)
tagsDown1_on_left_pad=[et[2] for et in entitiesDown1_on_left_pad[5:9]] 
#точки можно брать entities0_on_left_pad[6:10]
tagsDown0_on_left_pad=[et[2] for et in entities0_on_left_pad[6:10]] #индексы точек левой нижней панели
#создаем физические группы с одинаковым названием для левой верхней грани
gmsh.model.addPhysicalGroup(2, [d2[1][2]], -1, "DownLeft")
gmsh.model.addPhysicalGroup(1, tagsUp1_on_left_pad, -1, "DownLeft")
gmsh.model.addPhysicalGroup(0, tagsUp0_on_left_pad, -1, "DownLeft")



#Право верх
bb_UpRight = gmsh.model.getBoundingBox(d2[3]...)
entitiesUp1_on_right_pad = gmsh.model.getEntitiesInBoundingBox(bb_UpRight..., 1) #а справа все нормально почему то...
tagsUp1_on_right_pad = [et[2] for et in entitiesUp1_on_right_pad[1:5]]
#Точки право верх
entitiesUp0_on_right_pad = gmsh.model.getEntitiesInBoundingBox(bb_UpRight..., 0) 
tagsUp0_on_right_pad = [et[2] for et in entitiesUp0_on_right_pad[1:5]]
#создаем физические группы с одинаковым названием для правой верхней грани
gmsh.model.addPhysicalGroup(2, [d2[3][2]], -1, "UpRight")
gmsh.model.addPhysicalGroup(1, tagsUp1_on_right_pad, -1, "UpRight")
gmsh.model.addPhysicalGroup(0, tagsUp0_on_right_pad, -1, "UpRight")

#Право низ
#Новый бокс можно не создавать возьмем bb_UpRight[6:10]
tagsDown1_on_right_pad = [et[2] for et in entitiesUp1_on_right_pad[6:10]]
#Точки право низ
tagsDown0_on_right_pad = [et[2] for et in entitiesUp0_on_right_pad[6:10]]
#создаем физические группы с одинаковым названием для правой нижней грани
gmsh.model.addPhysicalGroup(2, [d2[13][2]], -1, "DownRight")
gmsh.model.addPhysicalGroup(1, tagsUp1_on_right_pad, -1, "DownRight")
gmsh.model.addPhysicalGroup(0, tagsUp0_on_right_pad, -1, "DownRight")




d3 = gmsh.model.getEntities(3)
gmsh.model.addPhysicalGroup(3, [d3[1][2]], -1, "Down")
gmsh.model.addPhysicalGroup(3, [d3[2][2]], -1, "Up")
gmsh.model.occ.synchronize
gmsh.model.occ.synchronize()

# gmsh.option.setNumber("Mesh.MeshSizeMin", 2.)
# gmsh.option.setNumber("Mesh.MeshSizeMax", 3.)
gmsh.model.mesh.generate(3)

msh_file = "model_balka.msh" |> tcf;
gmsh.write(msh_file)

# gmsh.fltk.run()

gmsh.finalize()




model = GmshDiscreteModel(msh_file);

vtk_file = "model_balka" |> tcf;
writevtk(model, vtk_file);




function get_nodes_by_tag(model::DiscreteModel, tag::String)::Vector{Int64}

    
    labels = get_face_labeling(model)
    tag_id = get_tag_from_name(labels, tag);
    
    dim = 0 # we will look for nodes which has dimension of 0
    
    dface_to_entity = get_face_entity(labels, dim) 
    
    dface_to_isontag = BitVector(undef, num_faces(labels,dim))

    tag_entities = get_tag_entities(labels, tag)

    for i in eachindex(dface_to_entity)
        buf = false
        for entity in tag_entities
            buf += dface_to_entity[i] == entity
        end
        dface_to_isontag[i] = buf
    end

    return findall(dface_to_isontag)

end


function get_nodes_in_bounding_box(model, bounding_box)
    node_coordinates = get_node_coordinates(model)
    nodes_in_bbox = BitVector(undef, length(node_coordinates))

    for i in eachindex(nodes_in_bbox)
        nodes_in_bbox[i] = 
            node_coordinates[i][1] >= bounding_box[1] && 
            node_coordinates[i][1] <= bounding_box[4] &&
            
            node_coordinates[i][2] >= bounding_box[2] && 
            node_coordinates[i][2] <= bounding_box[5] &&
            
            node_coordinates[i][3] >= bounding_box[3] && 
            node_coordinates[i][3] <= bounding_box[6];
    end
    return findall(nodes_in_bbox)
end

node_coordinates = get_node_coordinates(model)



# номера узлов граней. 
nodes_UpLeft = get_nodes_by_tag(model, "UpLeft")
nodes_UpRight = get_nodes_by_tag(model, "UpRight") #функция передает одинаковые узлы для верхнего Up и нижнего Down тега с каждой стороны.
#Причем передает в обоих случаях верхние узлы
nodes_DownLeft = setdiff(get_nodes_in_bounding_box(model, bb_UpLeft), nodes_UpLeft)
nodes_DownRight = setdiff(get_nodes_in_bounding_box(model, bb_UpRight), nodes_UpRight)

# coordinates_UpLeft = Vector{VectorValue{3, Float64}}()
coordinates_UpLeft = Vector{VectorValue{3, Any}}()
coordinates_DownLeft = Vector{VectorValue{3, Any}}()
j_UpLeft = Vector{Float64}()
j_DownLeft = Vector{Float64}()
for i in 1:length(nodes_UpLeft)
    push!(coordinates_UpLeft, node_coordinates[nodes_UpLeft[i]])
    push!(j_UpLeft, norm(coordinates_UpLeft[i]))

    push!(coordinates_DownLeft, node_coordinates[nodes_DownLeft[i]])
    push!(j_DownLeft, norm(coordinates_DownLeft[i]))
end
sorted_indices_j_UpLeft = sortperm(j_UpLeft)
sorted_j_UpLeft = sort(j_UpLeft)

sorted_indices_j_DownLeft = sortperm(j_DownLeft)
sorted_j_DownLeft = sort(j_DownLeft)


sorted_nodes_UpLeft = nodes_UpLeft[sorted_indices_j_UpLeft]
sorted_nodes_DownLeft = nodes_DownLeft[sorted_indices_j_DownLeft]
#sorted_nodes_DownLeft и sorted_nodes_UpLeft содержат глобальные индексы узлов верхней левой и нижней поверхности с одинаковыми координатами

#Также для правой части

coordinates_UpRight = Vector{VectorValue{3, Any}}()
coordinates_DownRight = Vector{VectorValue{3, Any}}()
j_UpRight = Vector{Float64}()
j_DownRight = Vector{Float64}()
for i in 1:length(nodes_UpRight)
    push!(coordinates_UpRight, node_coordinates[nodes_UpRight[i]])
    push!(j_UpRight, norm(coordinates_UpRight[i]))

    push!(coordinates_DownRight, node_coordinates[nodes_DownRight[i]])
    push!(j_DownRight, norm(coordinates_DownRight[i]))
end
sorted_indices_j_UpRight = sortperm(j_UpRight)
sorted_j_UpRight = sort(j_UpRight)

sorted_indices_j_DownRight = sortperm(j_DownRight)
sorted_j_DownRight = sort(j_DownRight)

sorted_nodes_UpRight = nodes_UpRight[sorted_indices_j_UpRight]
sorted_nodes_DownRight = nodes_DownRight[sorted_indices_j_DownRight]