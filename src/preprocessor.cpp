#include "preprocessor.h"
#include "global_variables.h"
#include <fstream>
#include "tinyxml2.h"
#include <sstream>
#include "initial_conditions.h"
#include "quad_bcs_plus.h"

using namespace tinyxml2;

preprocessor::preprocessor()
{
    //ctor

}

preprocessor::~preprocessor()
{
    //dtor
}

void preprocessor::initialise_program_variables(char* xmltest, global_variables &globals, domain_geometry &geometry,
                                                initial_conditions &initial_conds, quad_bcs_plus &bcs ) {

    XMLDocument xmlDoc;

    XMLError eResult = xmlDoc.LoadFile(xmltest);


    parse_global_variables(xmlDoc,  globals);
    parse_geometry_variables(xmlDoc,  geometry);
    parse_initial_conditions(xmlDoc, initial_conds);
    parse_boundary_conditions(xmlDoc, bcs);

    mach_number_factor(globals,bcs,initial_conds,geometry);
    geometry.scale_geometries(globals.scale);
    globals.initialise(geometry,initial_conds);



}
void preprocessor::mach_number_factor( global_variables &globals,quad_bcs_plus &bcs,
        initial_conditions &initials,domain_geometry &geometry ){


    double factor = globals.max_velocity/globals.scale ;
    //double factor = globals.max_velocity;

    bcs.e_u = bcs.e_u * factor;
    bcs.w_u = bcs.w_u * factor;
    bcs.n_u = bcs.n_u * factor;
    bcs.s_u = bcs.s_u *factor;

    bcs.e_v = bcs.e_v * factor;
    bcs.w_v = bcs.w_v * factor;
    bcs.n_v = bcs.n_v * factor;
    bcs.s_v = bcs.s_v *factor;

    initials.vel_gradient.x = initials.vel_gradient.x * factor ;
    initials.vel_gradient.y = initials.vel_gradient.y * factor;
    initials.vel_gradient.z = initials.vel_gradient.z * factor;
    initials.vel_origin_mag.x = initials.vel_origin_mag.x * factor;
    initials.vel_origin_mag.y = initials.vel_origin_mag.y * factor;
    initials.vel_origin_mag.z = initials.vel_origin_mag.z * factor;

    initials.velocity.x = initials.velocity.x * factor;
    initials.velocity.y = initials.velocity.y * factor;
    initials.velocity.z = initials.velocity.z * factor;

}
void preprocessor::parse_boundary_conditions(XMLDocument &xmlDoc, quad_bcs_plus &bcs){

    const char* parent = "boundary_conditions";



    /// WEST BCS
    bcs.w_rho = get_xml_double(parent,"west","rho",xmlDoc);
    bcs.w_u = get_xml_double(parent,"west","u",xmlDoc);
    bcs.w_v = get_xml_double(parent,"west","v",xmlDoc);
    std::string temp;

    temp = get_xml_text(parent,"west","vel_type",xmlDoc);
    if (temp.compare("dirichlet") == 0){
        bcs.w_type_vel = 1;
    }else if ( temp.compare("neumann") == 0){
        bcs.w_type_vel = 2;
    }else if ( temp.compare("periodic") == 0){
        bcs.w_type_vel = 3;
    }else if ( temp.compare("parabolic-N") == 0){
        bcs.w_type_vel = 4;
    }else if ( temp.compare("parabolic-W") == 0){
        bcs.w_type_vel = 5;
    }



    temp = get_xml_text(parent,"west","rho_type",xmlDoc);
    if (temp.compare("dirichlet") == 0){
        bcs.w_type_rho = 1;
    }else if ( temp.compare("neumann") == 0){
        bcs.w_type_rho = 2;
    }else if ( temp.compare("periodic") == 0){
        bcs.w_type_rho = 3;
    }

    /// east BCS
    bcs.e_rho = get_xml_double(parent,"east","rho",xmlDoc);
    bcs.e_u = get_xml_double(parent,"east","u",xmlDoc);
    bcs.e_v = get_xml_double(parent,"east","v",xmlDoc);


    temp = get_xml_text(parent,"east","vel_type",xmlDoc);
    if (temp.compare("dirichlet") == 0){
        bcs.e_type_vel = 1;
    }else if ( temp.compare("neumann") == 0){
        bcs.e_type_vel = 2;
    }else if ( temp.compare("periodic") == 0){
        bcs.e_type_vel = 3;
    }else if ( temp.compare("parabolic-N") == 0){
        bcs.w_type_vel = 4;
    }else if ( temp.compare("parabolic-W") == 0){
        bcs.w_type_vel = 5;
    }



    temp = get_xml_text(parent,"east","rho_type",xmlDoc);
    if (temp.compare("dirichlet") == 0){
        bcs.e_type_rho = 1;
    }else if ( temp.compare("neumann") == 0){
        bcs.e_type_rho = 2;
    }else if ( temp.compare("periodic") == 0){
        bcs.e_type_rho = 3;
    }

    /// north BCS
    bcs.n_rho = get_xml_double(parent,"north","rho",xmlDoc);
    bcs.n_u = get_xml_double(parent,"north","u",xmlDoc);
    bcs.n_v = get_xml_double(parent,"north","v",xmlDoc);


    temp = get_xml_text(parent,"north","vel_type",xmlDoc);
    if (temp.compare("dirichlet") == 0){
        bcs.n_type_vel = 1;
    }else if ( temp.compare("neumann") == 0){
        bcs.n_type_vel = 2;
    }else if ( temp.compare("periodic") == 0){
        bcs.n_type_vel = 3;
    }else if ( temp.compare("parabolic-N") == 0){
        bcs.w_type_vel = 4;
    }else if ( temp.compare("parabolic-W") == 0){
        bcs.w_type_vel = 5;
    }



    temp = get_xml_text(parent,"north","rho_type",xmlDoc);
    if (temp.compare("dirichlet") == 0){
        bcs.n_type_rho = 1;
    }else if ( temp.compare("neumann") == 0){
        bcs.n_type_rho = 2;
    }else if ( temp.compare("periodic") == 0){
        bcs.n_type_rho = 3;
    }

    /// south BCS
    bcs.s_rho = get_xml_double(parent,"south","rho",xmlDoc);
    bcs.s_u = get_xml_double(parent,"south","u",xmlDoc);
    bcs.s_v = get_xml_double(parent,"south","v",xmlDoc);


    temp = get_xml_text(parent,"south","vel_type",xmlDoc);
    if (temp.compare("dirichlet") == 0){
        bcs.s_type_vel = 1;
    }else if ( temp.compare("neumann") == 0){
        bcs.s_type_vel = 2;
    }else if ( temp.compare("periodic") == 0){
        bcs.s_type_vel = 3;
    }else if ( temp.compare("parabolic-N") == 0){
        bcs.w_type_vel = 4;
    }else if ( temp.compare("parabolic-W") == 0){
        bcs.w_type_vel = 5;
    }



    temp = get_xml_text(parent,"south","rho_type",xmlDoc);
    if (temp.compare("dirichlet") == 0){
        bcs.s_type_rho = 1;
    }else if ( temp.compare("neumann") == 0){
        bcs.s_type_rho = 2;
    }else if ( temp.compare("periodic") == 0){
        bcs.s_type_rho = 3;
    }

    // add in w velocity  afterwards for 3d
}



void preprocessor::parse_initial_conditions(XMLDocument &xmlDoc, initial_conditions &initials){

    const char* parent = "initial_conditions";
    initials.average_rho = get_xml_double(parent, "average_rho", xmlDoc);
    initials.pressure_gradient = get_xml_double(parent, "pressure_gradient", xmlDoc);
    initials.rho_gradient.x = get_xml_double(parent,"rho","gradient","x",xmlDoc);
    initials.rho_gradient.y = get_xml_double(parent,"rho","gradient","y",xmlDoc);
    initials.rho_gradient.z = get_xml_double(parent,"rho","gradient","z",xmlDoc);
    initials.rho_origin_mag.x = get_xml_double(parent,"rho","origin_magnitude","x",xmlDoc);
    initials.rho_origin_mag.y = get_xml_double(parent,"rho","origin_magnitude","y",xmlDoc);
    initials.rho_origin_mag.z = get_xml_double(parent,"rho","origin_magnitude","z",xmlDoc);
    initials.origin_loc.x = get_xml_double(parent,"origin_location","x",xmlDoc);
    initials.origin_loc.y = get_xml_double(parent,"origin_location","y",xmlDoc);
    initials.origin_loc.z = get_xml_double(parent,"origin_location","z",xmlDoc);
    initials.vel_gradient.x = get_xml_double(parent,"vel","gradient","x",xmlDoc);
    initials.vel_gradient.y = get_xml_double(parent,"vel","gradient","y",xmlDoc);
    initials.vel_gradient.z = get_xml_double(parent,"vel","gradient","z",xmlDoc);
    initials.vel_origin_mag.x = get_xml_double(parent,"vel","origin_magnitude","x",xmlDoc);
    initials.vel_origin_mag.y = get_xml_double(parent,"vel","origin_magnitude","y",xmlDoc);
    initials.vel_origin_mag.z = get_xml_double(parent,"vel","origin_magnitude","z",xmlDoc);

    initials.velocity.x = get_xml_double(parent, "u", xmlDoc);
    initials.velocity.y = get_xml_double(parent, "v", xmlDoc);
    initials.velocity.z = get_xml_double(parent, "w", xmlDoc);

}
void preprocessor::parse_geometry_variables(XMLDocument &xmlDoc, domain_geometry &geometry){

    const char* parent = "geometry";
    geometry.X = get_xml_double(parent, "x", xmlDoc);
    geometry.Y= get_xml_double(parent, "y",xmlDoc);
    geometry.dx = get_xml_double(parent, "dx", xmlDoc);
    geometry.dy = get_xml_double(parent, "dy", xmlDoc);
    geometry.dt= get_xml_double(parent, "streaming_dt", xmlDoc);

    geometry.initialise();


}
void preprocessor::parse_global_variables(XMLDocument &xmlDoc, global_variables &globals){

    const char* parent = "global_variables";
    globals.tolerance = get_xml_double(parent, "tolerance", xmlDoc);
    globals.pre_conditioned_gamma = get_xml_double(parent, "pre_condition_gamma",xmlDoc);
    globals.simulation_length = get_xml_double(parent, "simulation_length", xmlDoc);
    globals.time_marching_step = get_xml_double(parent, "time_marching_step", xmlDoc);
    globals.reynolds_number = get_xml_double(parent, "reynolds_no", xmlDoc);
    globals.max_velocity = get_xml_double(parent, "mach_no", xmlDoc);
    globals.simulation_name = get_xml_text(parent,"simulation_name",xmlDoc);
    globals.output_file_dir = get_xml_text(parent, "output_directory", xmlDoc);
    globals.max_mg_levels = get_xml_double(parent,"max_multi_grid_levels",xmlDoc);
    globals.fmg_levels =  get_xml_double(parent,"FMG_levels",xmlDoc);
    globals.fmg_tolerance = get_xml_double(parent,"FMG_tolerance",xmlDoc);
    globals.arti_disp_kappa_2 = get_xml_double(parent,"dissipation_kappa_2",xmlDoc);
    globals.arti_disp_kappa_4 = get_xml_double(parent,"dissipation_kappa_4",xmlDoc);
    globals.martinelli = get_xml_double(parent,"martinelli_exponent",xmlDoc);
    globals.testcase= get_xml_double(parent, "testcase", xmlDoc);
    globals.mesh_type= get_xml_double(parent, "mesh_type", xmlDoc);
    globals.scale= get_xml_double(parent, "scale", xmlDoc);
    globals.womersley_no= get_xml_double(parent, "womersley_no", xmlDoc);
}

double preprocessor::get_xml_double(const char* parent, const char* child, XMLDocument &doc){

    XMLError eResult;
    XMLNode * pRoot = doc.FirstChild();

    if (pRoot == nullptr) {

        eResult = XML_ERROR_FILE_READ_ERROR;
    }

    XMLElement * parent_element = pRoot->FirstChildElement(parent);
    if (parent_element == nullptr) eResult =  XML_ERROR_PARSING_ELEMENT;
    XMLElement* child_element = parent_element -> FirstChildElement( child );
    if (child_element == nullptr) eResult =  XML_ERROR_PARSING_ELEMENT;

    double output;

    eResult = child_element->QueryDoubleText(&output);

    return output;
}
double preprocessor::get_xml_double(const char* parent, const char* child, const char* child2,
                                    const char* child3, XMLDocument &doc){

    XMLError eResult;
    XMLNode * pRoot = doc.FirstChild();

    if (pRoot == nullptr) {

        eResult = XML_ERROR_FILE_READ_ERROR;
    }

    XMLElement * parent_element = pRoot->FirstChildElement(parent);
    if (parent_element == nullptr) eResult =  XML_ERROR_PARSING_ELEMENT;
    XMLElement* child_element = parent_element -> FirstChildElement( child );
    if (child_element == nullptr) eResult =  XML_ERROR_PARSING_ELEMENT;
    XMLElement* child2_element = child_element -> FirstChildElement( child2 );
    if (child2_element == nullptr) eResult =  XML_ERROR_PARSING_ELEMENT;
    XMLElement* child3_element = child2_element -> FirstChildElement( child3 );
    if (child3_element == nullptr) eResult =  XML_ERROR_PARSING_ELEMENT;
    double output;

    eResult = child3_element->QueryDoubleText(&output);

    return output;
}

double preprocessor::get_xml_double(const char* parent, const char* child, const char* child2,
                                     XMLDocument &doc){

    XMLError eResult;
    XMLNode * pRoot = doc.FirstChild();

    if (pRoot == nullptr) {

        eResult = XML_ERROR_FILE_READ_ERROR;
    }

    XMLElement * parent_element = pRoot->FirstChildElement(parent);
    if (parent_element == nullptr) eResult =  XML_ERROR_PARSING_ELEMENT;
    XMLElement* child_element = parent_element -> FirstChildElement( child );
    if (child_element == nullptr) eResult =  XML_ERROR_PARSING_ELEMENT;
    XMLElement* child2_element = child_element -> FirstChildElement( child2 );
    if (child2_element == nullptr) eResult =  XML_ERROR_PARSING_ELEMENT;

    double output;

    eResult = child2_element->QueryDoubleText(&output);

    return output;
}
const char * preprocessor::get_xml_text(const char* parent, const char* child,XMLDocument &doc){

    XMLError eResult;
    XMLNode * pRoot = doc.FirstChild();

    if (pRoot == nullptr) {

        eResult = XML_ERROR_FILE_READ_ERROR;
    }


    XMLElement * parent_element = pRoot->FirstChildElement(parent);
    if (parent_element == nullptr) eResult =  XML_ERROR_PARSING_ELEMENT;
    XMLElement* child_element = parent_element -> FirstChildElement( child );
    if (child_element == nullptr) eResult =  XML_ERROR_PARSING_ELEMENT;

    const char * output;

    output = child_element->GetText();

    return output;
}
const char * preprocessor::get_xml_text(const char* parent, const char* child, const char* child2
                                        , XMLDocument &doc){

    XMLError eResult;
    XMLNode * pRoot = doc.FirstChild();

    if (pRoot == nullptr) {

        eResult = XML_ERROR_FILE_READ_ERROR;
    }


    XMLElement * parent_element = pRoot->FirstChildElement(parent);
    if (parent_element == nullptr) eResult =  XML_ERROR_PARSING_ELEMENT;
    XMLElement* child_element = parent_element -> FirstChildElement( child );
    if (child_element == nullptr) eResult =  XML_ERROR_PARSING_ELEMENT;
    XMLElement* child2_element = child_element -> FirstChildElement( child2 );
    if (child2_element == nullptr) eResult =  XML_ERROR_PARSING_ELEMENT;

    const char * output;

    output = child2_element->GetText();

    return output;
}
