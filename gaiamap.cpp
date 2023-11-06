/***

	gaiamap.cpp

	Use the 1 square degree data files to produce a custom star chart based on Gaia data.

	Copyright (c) 2022, Chris Street.

***/
#define CODEHAPPY_NATIVE
#include <libcodehappy.h>

#define DATA_DIR	"degrees/"

/* TrueType fonts that we use */
static Font swansea(font_swansea);
static Font oregon(font_oregon);
static Font italic(font_oregon_italic);

/* Options for the map render. */
struct RenderOptions {
	RenderOptions();
	// do we write right ascension as degrees/hours/minutes or hours/minutes/seconds?
	bool ra_deg;
	// the limiting magnitude (only draw sources brighter than this value)
	double mag_limit;
	// the color for designations
	RGBColor desig_color;
	// the color for coordinates
	RGBColor coord_color;
	// the color for the background
	RGBColor bg_color;
	// the color for the coordinate grid
	RGBColor grid_color;
	// do we dim sources that are far away?
	bool dim_far;
	// show decimal degrees for both RA and DEC (overrides ra_deg)
	bool decimal_degrees;
	// draw a key at the bottom of the chart
	bool draw_key;
	// draw distance in italics below the star's designation (intersections permitting)
	bool show_distance;
	// draw source designations on the chart
	bool show_desig;
	// distance limit (only draw sources closer than this value; <= 0. for no limit)
	double close_limit;
	// distance limit (only draw sources farther than this value; <= 0. for no limit)
	double far_limit;
	// always use the shortest designation for a source
	bool short_desig;
	// output the cone file?
	bool out_cone;
	// render the map?
	bool no_map;
	// use smaller stars?
	bool small_stars;
	// suppress UCAC4 designations
	bool no_ucac;
	// show SIMBAD "Cl*", "Cl", "V*", "EM*", "**", "NAME", etc. designations as-is
	bool simbad_asis;
	// use negative colors (white background, black stars) -- good for printing
	bool neg_printable;
	// use a strictly rectangonal projection (at high declinations, this will cause things to look 'squashed'.)
	bool rect_proj;
	// RA should not be adjusted for better aspect ratio: change the width/height instead.
	bool fix_ra;
	// turn off TIC (TESS Input Catalog) designations (we could match a better designation)
	bool use_tic;
	// a declination to use for cos(DEC) adjustment (use it to tile multiple non-rectangular projections.)
	double use_dec;
	// draw proper motion arrows?
	bool draw_pm;
	// the number of years of proper motion
	double pm_years;
	// do we draw a coordinate grid?
	bool draw_grid;
	// downsample the final image to this size (for photograph-style renders.)
	u32 downsample_width;
};

RenderOptions::RenderOptions() {
	ra_deg = false;
	mag_limit = 999999.;
	desig_color = C_YELLOW;
	coord_color = MAKE_RGB(0, 0xff, 0);
	bg_color = MAKE_RGB(4, 46, 96);
	dim_far = false;
	decimal_degrees = false;
	draw_key = false;
	show_distance = false;
	show_desig = true;
	close_limit = -1.0;
	far_limit = -1.0;
	short_desig = false;
	out_cone = true;
	no_map = false;
	small_stars = false;
	no_ucac = false;
	simbad_asis = false;
	neg_printable = false;
	grid_color = RGB_NO_CHECK(0xa0, 0xa0, 0xa0);
	rect_proj = false;
	fix_ra = false;
	use_tic = true;
	use_dec = 999.0;
	draw_pm = false;
	pm_years = 1000.0;
	draw_grid = true;
	downsample_width = 0;
}

static RenderOptions ro;

void bv2rgb(double &r,double &g,double &b,double bv)    // RGB <0,1> <- BV <-0.4,+2.0> [-]
    {
    double t;  r=0.0; g=0.0; b=0.0; if (bv<-0.4) bv=-0.4; if (bv> 2.0) bv= 2.0;
         if ((bv>=-0.40)&&(bv<0.00)) { t=(bv+0.40)/(0.00+0.40); r=0.61+(0.11*t)+(0.1*t*t); }
    else if ((bv>= 0.00)&&(bv<0.40)) { t=(bv-0.00)/(0.40-0.00); r=0.83+(0.17*t)          ; }
    else if ((bv>= 0.40)&&(bv<2.10)) { t=(bv-0.40)/(2.10-0.40); r=1.00                   ; }
         if ((bv>=-0.40)&&(bv<0.00)) { t=(bv+0.40)/(0.00+0.40); g=0.70+(0.07*t)+(0.1*t*t); }
    else if ((bv>= 0.00)&&(bv<0.40)) { t=(bv-0.00)/(0.40-0.00); g=0.87+(0.11*t)          ; }
    else if ((bv>= 0.40)&&(bv<1.60)) { t=(bv-0.40)/(1.60-0.40); g=0.98-(0.16*t)          ; }
    else if ((bv>= 1.60)&&(bv<2.00)) { t=(bv-1.60)/(2.00-1.60); g=0.82         -(0.5*t*t); }
         if ((bv>=-0.40)&&(bv<0.40)) { t=(bv+0.40)/(0.40+0.40); b=1.00                   ; }
    else if ((bv>= 0.40)&&(bv<1.50)) { t=(bv-0.40)/(1.50-0.40); b=1.00-(0.47*t)+(0.1*t*t); }
    else if ((bv>= 1.50)&&(bv<1.94)) { t=(bv-1.50)/(1.94-1.50); b=0.63         -(0.6*t*t); }
    }

/* From B and V magnitudes, return the RGBColor representing the blackbody. */
RGBColor bv_to_rgbcolor(double b, double v) {
	double rr, gg, bb;
	if (v >= 20. || b >= 20.)
		b = v;
	bv2rgb(rr, gg, bb, b - v);
	u32 r, g, bl;
	r = (u32)floor(255. * rr + 0.5);
	g = (u32)floor(255. * gg + 0.5);
	bl = (u32)floor(255. * bb + 0.5);
	return RGB_NO_CHECK(r, g, bl);
}


/* Given a magnitude, return the radius of the circle representing the star. */
u32 radius_for_magnitude(double mag) {
	const u32 rad_from_mag[] = {
		40, 36, 32, 28, 24, 22, 20, 18, 16, 14, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1
	//      -1   0   1   2   3   4   5   6   7   8   9  10  11 12 13 14 15 16 17 18 19 20
	};
	const u32 rad_from_mag_smol[] = {
		20, 18, 16, 14, 12, 10,  8,  7,  6,  5,  4,  3,  2, 1, 1, 1, 1, 1, 1, 1, 1, 1
	//      -1   0   1   2   3   4   5   6   7   8   9  10  11 12 13 14 15 16 17 18 19 20
	};
	int idx = floor(mag) + 1;
	if (idx < 0)
		idx = 0;
	if (idx > 21)
		idx = 21;
	return ro.small_stars ? rad_from_mag_smol[idx] : rad_from_mag[idx];
}

u32 font_size_for_height(Font& font, u32 height) {
	SBitmap* test;
	u32 ret;

	test = font.render_cstr("test string", 200, false, nullptr);
	ret = (200 * height + (test->height() / 2)) / test->height();
	delete test;

	return ret;
}

u32 font_size_for_height_cached(Font& font, u32 height) {
	u32 cache[400] = { 0 };
	if (height >= 400)
		return font_size_for_height(font, height);
	if (cache[height] == 0)
		cache[height] = font_size_for_height(font, height);
	return cache[height];
}

struct MapData {
	std::string desig[6];
	std::string category;
	std::string spectral_type;
	u64 s_id;
	double ra;
	double dec;
	float mag_g;
	float mag_bp;
	float mag_rp;
	float dist;
	float pmra;
	float pmdec;
	float radvel;
	double angsize[2];
	double l;
	double b;
	double MagnitudeWeighted(void) const;
	bool IsGalaxy(void) const;
	u32 DataScore(void) const;
	void ClearDesig(void);
};

u32 MapData::DataScore(void) const {
	u32 ret = 0;
	if (mag_g < 20.)
		ret += 2;
	if (mag_bp < 20.)
		++ret;
	if (mag_rp < 20.)
		++ret;
	if (dist > 0.)
		ret += 5;
	if (pmra != 0. && pmdec != 0.)
		ret += 5;
	if (radvel != 0.)
		ret += 2;
	if (angsize[0] != 0.)
		++ret;
	return ret;
}

void MapData::ClearDesig(void) {
	for (int e = 0; e < 6; ++e)
		desig[e].clear();
	category.clear();
	spectral_type.clear();
}

double MapData::MagnitudeWeighted(void) const {
	double magw = 1.0;
	u32 cmag = 0;
	if (mag_bp <= 20.) {
		magw *= exp(mag_bp);
		cmag++;
	}
	if (mag_rp <= 20.) {
		magw *= exp(mag_rp);
		cmag++;
	}
	if (mag_g <= 20.) {
		magw *= exp(mag_g);
		cmag++;
	}
	if (0 == cmag)
		return 20.;
	magw = pow(magw, 1.0 / (double)cmag);
	magw = log(magw);
	return magw;
}

bool MapData::IsGalaxy(void) const {
	const char* glx_cats[] = {
	"Active Galaxy Nucleus", "Blue compact Galaxy", "Brightest galaxy in a Cluster (BCG)", "Cluster of Galaxies",
	"Emission-line galaxy", "Galaxy", "Galaxy in Cluster of Galaxies", "Galaxy in Group of Galaxies", "Galaxy in Pair of Galaxies",
	"HII Galaxy", "Interacting Galaxies", "LINER-type Active Galaxy Nucleus", "Low Surface Brightness Galaxy",
	"Pair of Galaxies", "Possible Active Galaxy Nucleus", "Possible Quasar", "Quasar", "Radio Galaxy", "Seyfert 1 Galaxy",
	"Seyfert 2 Galaxy", "Seyfert Galaxy", "Starburst Galaxy"
	};
	for (auto n : glx_cats)
		if (category == n)
			return true;
	if (!desig[0].empty() && !strncmp(desig[0].c_str(), "LEDA ", 5))
		return true;
	if (!desig[0].empty() && !strncmp(desig[0].c_str(), "Gaia-UNWISE ", 12))
		return true;
	if (!category.empty() && strstr(category.c_str(), "Galaxy") != nullptr)
		return true;
	if (!category.empty() && strstr(category.c_str(), "galaxy") != nullptr)
		return true;
	return false;
}

bool mapdata_comp(const MapData& md1, const MapData& md2) {
	double mag1, mag2;
	mag1 = std::min(md1.mag_g, md1.mag_bp);
	mag2 = std::min(md2.mag_g, md2.mag_bp);
	return mag1 < mag2;
}

bool mapdata_compa(const MapData& md1, const MapData& md2) {
	double mag1, mag2;
	mag1 = std::min(md1.mag_g, md1.mag_bp);
	mag2 = std::min(md2.mag_g, md2.mag_bp);
	return mag1 > mag2;
}

bool mapdata_compn(const MapData& md1, const MapData& md2) {
	return md1.desig[0] < md2.desig[0];
}

struct Rect {
	Rect(int x1_, int x2_, int y1_, int y2_);
	bool in(int x, int y) const;
	bool isect(int x1_, int x2_, int y1_, int y2_) const;
	int x1;
	int x2;
	int y1;
	int y2;
};

Rect::Rect(int x1_, int x2_, int y1_, int y2_) {
	SORT2(x1_, x2_, int);
	SORT2(y1_, y2_, int);
	x1 = x1_;
	x2 = x2_;
	y1 = y1_;
	y2 = y2_;
}

bool Rect::in(int x, int y) const {
	return (x >= x1 && x <= x2 && y >= y1 && y <= y2);
}

bool Rect::isect(int x1_, int x2_, int y1_, int y2_) const {
	SORT2(x1_, x2_, int);
	SORT2(y1_, y2_, int);
	x1_ = std::max(x1_, x1);
	x2_ = std::min(x2_, x2);
	y1_ = std::max(y1_, y1);
	y2_ = std::min(y2_, y2);
	return (x1_ <= x2_ && y1_ <= y2_);
}

class Intersection {
public:
	Intersection(u32 width, u32 height);
	void add_rect(int x1, int x2, int y1, int y2);
	bool test_rect(int x1, int x2, int y1, int y2) const;
	void clear();

private:
	std::vector<Rect> rects;
};

Intersection::Intersection(u32 width, u32 height) {
	clear();
}

void Intersection::clear(void) {
	rects.clear();
}

void Intersection::add_rect(int x1, int x2, int y1, int y2) {
	rects.push_back(Rect(x1, x2, y1, y2));
}

bool Intersection::test_rect(int x1, int x2, int y1, int y2) const {
	SORT2(x1, x2, int);
	SORT2(y1, y2, int);
	for (const auto& r : rects) {
		if (r.isect(x1, x2, y1, y2))
			return false;
	}
	return true;
}

RGBColor dim_color_by_distance(RGBColor clr, float dist) {
	float rr = (float)RGB_RED(clr),
	      gg = (float)RGB_GREEN(clr),
	      bb = (float)RGB_BLUE(clr);
	float intsy;
	if (dist <= 0.)
		dist = 5000.;
	if (dist <= 500.)
		intsy = 1.;
	else if (dist >= 5000.)
		intsy = 0.6;
	else
		intsy = 1.0 - 0.4 * ((dist - 500.) / 4500.);
	rr *= intsy;
	gg *= intsy;
	bb *= intsy;
	rr = floor(rr + 0.5);
	gg = floor(gg + 0.5);
	bb = floor(bb + 0.5);
	return MAKE_RGB((int)rr, (int)gg, (int)bb);
}

class Map {
public:
	Map(double ra1, double ra2, double dec1, double dec2, u32 width, u32 height);
	void render(void);
	void to_file(const char* fname);
	void render_zoom_in(void);

private:
	void read_regional_map_data(void);
	void read_map_data_file(const char* fname);
	void ra_dec_to_xy(double ra_in, double dec_in, int& x_out, int& y_out) const;
	void render_name(const std::string& name, const std::string& spectral_type, int x, int y, u32 rad, double dist);
	double step_value(double val1, double val2);
	void out_cone(void);
	void out_cluster_csv(void);

	double ra[2], dec[2];
	double cos_ra, rra[2];
	u32 w, h;
	std::vector<MapData> map_data;
	SBitmap* map;
	Intersection* isect;
};

double deg_to_rad(double deg) {
	deg *= 3.1415926535;
	deg /= 180.;
	return deg;
}

Map::Map(double ra1, double ra2, double dec1, double dec2, u32 width, u32 height) {
	SORT2(ra1, ra2, double);
	SORT2(dec1, dec2, double);
	ra[0] = ra1;
	ra[1] = ra2;
	dec[0] = dec1;
	dec[1] = dec2;
	if (ra[0] < 0.)
		ra[0] = 0.;
	if (ra[1] > 360.)
		ra[1] = 360.;
	if (dec[0] < -90.)
		dec[0] = -90.;
	if (dec[1] > 90.)
		dec[1] = 90.;
	w = width;
	h = height;
	if (ro.use_dec < 90. && ro.use_dec > -90.)
		cos_ra = cos(deg_to_rad(ro.use_dec));
	else
		cos_ra = cos(deg_to_rad((dec[0] + dec[1]) * 0.5));
	bool show_dimensions = (0 == h || 0 == w);
	if (!ro.rect_proj && ro.fix_ra && h != 0 && w != 0) {
		std::cout << "*** Warning: you should constrain only one of width or height when -fixra option is on.\n";
		std::cout << "[Switching to a rectagonal projection to preserve both RA range and specified dimensions.]\n";
		ro.fix_ra = false;
		ro.rect_proj = true;
	}
	if (0 == h) {
		if (0 == w)
			w = 4000;
		double ar = (dec[1] - dec[0]) / (ra[1] - ra[0]);
		if (!ro.rect_proj && ro.fix_ra) {
			h = (u32)floor(ar * w / cos_ra + 0.5);
		} else {
			h = (u32)floor(ar * w + 0.5);
		}
	} else if (0 == w) {
		double ar = (ra[1] - ra[0]) / (dec[1] - dec[0]);
		if (!ro.rect_proj && ro.fix_ra) {
			w = (u32)floor(ar * h * cos_ra + 0.5);
		} else {
			w = (u32)floor(ar * h + 0.5);
		}
	}
	if (show_dimensions)
		std::cout << "Chart dimensions: " << w << " x " << h << std::endl;
	if (ro.rect_proj || ro.fix_ra) {
		rra[0] = ra[0];
		rra[1] = ra[1];
	} else {
		double mra = (ra[0] + ra[1]) * 0.5;
		double wra = (ra[1] - ra[0]) * 0.5;
		wra /= cos_ra;
		rra[0] = mra - wra;
		rra[1] = mra + wra;
		if (rra[0] < 0.)
			rra[0] = 0.;
		if (rra[1] > 360.)
			rra[1] = 360.;
		std::cout << ">>> RA range adjusted to [" << rra[0] << "," << rra[1] << "] for more accurate aspect ratio.\n";
		std::cout << "[You can turn off this adjustment and use RA strictly as entered with the option -rect.]\n";
		std::cout << "[You can also force an adjustment in width or height instead of RA with -fixra.]\n";
	}
	if (ro.neg_printable)
		map = new SBitmap(w, h, BITMAP_GRAYSCALE);
	else
		map = new SBitmap(w, h);
	isect = new Intersection(w, h);
	read_regional_map_data();
	std::cout << ">>> " << map_data.size() << " total sources read from disk.\n";
}

void Map::read_regional_map_data(void) {
	int x1, x2, y1, y2;
	x1 = (int)floor(rra[0]);
	x2 = (int)floor(rra[1]);
	if (rra[1] - floor(rra[1]) == 0.0)
		--x2;
	if (x2 < x1)
		x2 = x1;
	y1 = (int)floor(dec[0] + 90.0);
	y2 = (int)floor(dec[1] + 90.0);
	if (dec[1] - floor(dec[1]) == 0.0)
		--y2;
	if (y2 < y1)
		y2 = y1;
	for (int e = x1; e <= x2; ++e)
		for (int f = y1; f <= y2; ++f) {
			char fname[64];
			sprintf(fname, DATA_DIR "%06d", e * 1000 + f);
			read_map_data_file(fname);
		}

	// Let's check for duplicated designations; this sometimes happens in the SIMBAD Vizier cross match dump.
	std::sort(map_data.begin(), map_data.end(), mapdata_compn);
	for (int e = 1; e < map_data.size(); ++e) {
		if (map_data[e - 1].desig[0] != map_data[e].desig[0])
			continue;
		// These two have the same designation... preserve the one that scores higher on Gaia data, or comes first on a tie.
		u32 s1 = map_data[e - 1].DataScore(), s2 = map_data[e].DataScore();
		if (s1 > s2) {
			map_data[e].ClearDesig();
		} else if (s2 > s1) {
			map_data[e - 1].ClearDesig();
		} else {
			map_data[e].ClearDesig();
		}
	}
}

void Map::read_map_data_file(const char* fname) {
	std::ifstream i;
	char line[1024];
	std::cout << "Reading map data from " << fname << "...\n";

	i.open(fname);
	forever {
		MapData md;
		char* w, * w2;

		i.getline(line, 1024);
		if (i.eof())
			break;

		md.s_id = strtoull(line, nullptr, 10);
		w = strchr(line, '|');
		if (is_null(w))
			break;
		++w;
		// Common designation
		w2 = strchr(w, '|');
		*w2 = 0;
		md.desig[0] = w;
		w = w2 + 1;

		// Other designations
		for (int e = 1; e < 6; ++e) {
			w2 = strchr(w, '|');
			*w2 = 0;
			md.desig[e] = w;
			w = w2 + 1;
		}

		// SIMBAD primary category
		w2 = strchr(w, '|');
		*w2 = 0;
		md.category = w;
		w = w2 + 1;

		// SIMBAD spectral type
		w2 = strchr(w, '|');
		*w2 = 0;
		md.spectral_type = w;
		w = w2 + 1;

		// RA/DEC
		md.ra = atof(w);
		w = strchr(w, '|') + 1;
		md.dec = atof(w);
		w = strchr(w, '|') + 1;

		// Magnitude (G/V-band)
		if (*w == '|')
			md.mag_g = 40.;
		else
			md.mag_g = atof(w);
		w = strchr(w, '|') + 1;
		// B-band
		if (*w == '|')
			md.mag_bp = 40.;
		else
			md.mag_bp = atof(w);
		w = strchr(w, '|') + 1;
		// R-band
		if (*w == '|')
			md.mag_rp = 40.;
		else
			md.mag_rp = atof(w);
		w = strchr(w, '|') + 1;

		// Check against the magnitude limit, if any.
		if (md.mag_g > ro.mag_limit && md.mag_bp > ro.mag_limit && md.mag_rp > ro.mag_limit)
			continue;

		// Distance (converted from parsecs to light years)
		if (*w == '|')
			md.dist = -1.;
		else
			md.dist = atof(w) * 3.261564;
		w = strchr(w, '|') + 1;

		// Check against the distance limits, if any.
		if (ro.close_limit > 0. && md.dist > ro.close_limit)
			continue;
		if (ro.far_limit > 0. && md.dist < ro.far_limit)
			continue;

		// Proper motion
		md.pmra = atof(w);
		w = strchr(w, '|') + 1;
		md.pmdec = atof(w);
		w = strchr(w, '|') + 1;

		// Angular size.
		md.angsize[0] = atof(w);
		w = strchr(w, '|') + 1;
		md.angsize[1] = atof(w);
		w = strchr(w, '|') + 1;

		// Galactic longitude and latitude.
		md.l = atof(w);
		w = strchr(w, '|') + 1;
		md.b = atof(w);
		w = strchr(w, '|') + 1;

		// Radial velocity.
		md.radvel = atof(w);

		// Store it.
		map_data.push_back(md);
	}
	i.close();
	i.clear();
}

void Map::render(void) {
//	bool labels = yes_no("Include star designations?");
	bool labels = ro.show_desig;
	if (ro.no_map)
		return;
	std::cout << "Rendering map...\n";
	map->clear(ro.bg_color);

	// First, plot the coordinate grid.
	int x, y;
	double step = step_value(rra[0], rra[1]);
	double coor = floor(rra[0] / step) * step;
	int w_ra, w_dec, mx = 9999999, my = 9999999, _unused;
	ra_dec_to_xy(rra[1] - step, 0, w_ra, _unused);
	forever {
		ra_dec_to_xy(coor, dec[1], x, y);
		if (x < 0)
			break;
		if (x < mx)
			mx = x;
		coor += step;
	}
	step = step_value(dec[0], dec[1]);
	ra_dec_to_xy(0, dec[1] - step, _unused, w_dec);
	coor = floor(dec[0] / step) * step;
	forever {
		ra_dec_to_xy(ra[1], coor, x, y);
		if (y < 0)
			break;
		if (y < my)
			my = y;
		coor += step;
	}

	step = step_value(rra[0], rra[1]);
	coor = floor(rra[0] / step) * step;
	if (ro.draw_grid)
	   forever {
		int deg, min, sec;
		double rem;
		char pos[32];
		ra_dec_to_xy(coor, dec[1], x, y);
		if (x < 0)
			break;
		for (y = 0; y < h; ++y) {
			if (y % 5 > 2)
				continue;
			map->put_pixel(x, y, ro.grid_color);
		}
		deg = (int)floor(std::abs(coor));
		if (coor < 0.)
			deg = -deg;
		rem = std::abs(coor) - floor(std::abs(coor));
		min = (int)floor(60.0 * rem);
		rem -= double(min) / 60.0;
		sec = (int)floor(3600.0 * rem + 0.5);
		if (60 == sec) {
			++min;
			sec = 0;
		}
		if (60 == min) {
			if (deg < 0)
				--deg;
			else
				++deg;
			min = 0;
		}
		sprintf(pos, "%d  %02d' %02d\"", deg, min, sec);
		int idx = 0;
		while (!isspace(pos[idx]))
			++idx;
		ustring ustr = cstr2ustr(pos);
		ustr[idx] = 0xB0;
		SBitmap* render = oregon.render_ustr(ustr, 18, false, nullptr);
		delete [] ustr;
		for (y = my + w_dec / 2; y < h; y += w_dec) {
			Font::blit(render, map, x + 4, y - render->height() / 2, ro.coord_color);
		}
		delete render;
		coor += step;
	   }
	step = step_value(dec[0], dec[1]);
	coor = floor(dec[0] / step) * step;
	if (ro.draw_grid)
	   forever {
		int deg, min, sec;
		double rem;
		char pos[32];
		ra_dec_to_xy(ra[1], coor, x, y);
		if (y < 0)
			break;
		for (x = 0; x < w; ++x) {
			if (x % 5 > 2)
				continue;
			map->put_pixel(x, y, ro.grid_color);
		}
		deg = (int)floor(std::abs(coor));
		if (coor < 0.)
			deg = -deg;
		rem = std::abs(coor) - floor(std::abs(coor));
		min = (int)floor(60.0 * rem);
		rem -= double(min) / 60.0;
		sec = (int)floor(3600.0 * rem + 0.5);
		if (60 == sec) {
			++min;
			sec = 0;
		}
		if (60 == min) {
			if (deg < 0)
				--deg;
			else
				++deg;
			min = 0;
		}
		sprintf(pos, "%d  %02d' %02d\"", deg, min, sec);
		int idx = 0;
		while (!isspace(pos[idx]))
			++idx;
		ustring ustr = cstr2ustr(pos);
		ustr[idx] = 0xB0;
		SBitmap* render = oregon.render_ustr(ustr, 18, false, nullptr);
		delete [] ustr;
		for (x = mx + w_ra / 2; x < w; x += w_ra) {
			Font::blit(render, map, x - render->width() / 2, y + 4, ro.coord_color);
		}
		delete render;
		coor += step;
	   }

	// Now, plot the stars, dimmest to brightest.
	std::sort(map_data.begin(), map_data.end(), mapdata_compa);
	for (const auto& md : map_data) {
		int x, y;
		ra_dec_to_xy(md.ra, md.dec, x, y);
		if (x < 0 || x >= w)
			continue;
		if (y < 0 || y >= h)
			continue;
		double mag_b, mag_g;
		mag_b = md.mag_bp;
		mag_g = md.mag_g;
		if (mag_b >= 30.)
			mag_b = mag_g;
		if (mag_g >= 30.) {
			mag_b = 20.;
			mag_g = 20.;
		}
		u32 rad = radius_for_magnitude(md.MagnitudeWeighted());
		RGBColor star_color = bv_to_rgbcolor(mag_b, mag_g);
		if (ro.neg_printable)
			star_color = C_BLACK;
		else if (ro.dim_far)
			star_color = dim_color_by_distance(star_color, md.dist);
		if (md.IsGalaxy()) {
			// plot galaxies with 'X'es
			if (rad < 2)
				rad = 2;
			map->line(x - rad, y - rad, x + rad, y + rad, star_color);
			map->line(x + rad, y - rad, x - rad, y + rad, star_color);
		} else {
			map->fillcircle(x, y, rad, star_color);
		}
	}

	// Render proper motion arrows.
	for (const auto& md : map_data) {
		int x1, y1, x2, y2;
		if (!ro.draw_pm)
			break;
		if (md.pmra == 0. && md.pmdec == 0.)
			continue;
		// 1000 years of proper motion
		ra_dec_to_xy(md.ra, md.dec, x1, y1);
		ra_dec_to_xy(md.ra + (md.pmra * ro.pm_years) / 3600000., md.dec + (md.pmdec * ro.pm_years) / 3600000., x2, y2);
		map->line(x1, y1, x2, y2, C_RED);
	}

	// Now, render the designations. We use an Intersection object to avoid a bunch of overlapping labels, and
	// render from brightest star to dimmest star to favor brighter star designations.
	if (labels) {
		std::sort(map_data.begin(), map_data.end(), mapdata_comp);
		for (const auto& md : map_data) {
			if (md.desig[0].empty())
				continue;
			const std::string* use = &md.desig[0];
			if (!strncmp(md.desig[0].c_str(), "Gaia ", 5)) {
				if (!md.desig[5].empty())
					use = &md.desig[5];
			}
			if (ro.short_desig && !md.desig[5].empty()) {
				use = &md.desig[5];
			}
			if (!ro.use_tic && !strncmp(use->c_str(), "TIC ", 4)) {
				continue;
			}
			int x, y;
			ra_dec_to_xy(md.ra, md.dec, x, y);
			u32 rad = radius_for_magnitude(md.MagnitudeWeighted());
			double d = md.dist;
			if (md.IsGalaxy())
				d = -1.;
			render_name(*use, md.spectral_type, x, y, rad, d);
		}
	}
}

void Map::render_zoom_in(void) {
	double ra_c, dec_c, ra_w, dec_h;
	u32 fn = 0;
	char fname[64];

	// Default draw options.
	ro.draw_grid = false;
	ro.show_desig = false;
	ro.draw_pm = false;
	ro.show_distance = false;
	ro.bg_color = C_BLACK;
	ro.out_cone = false;

	// Center and radius of zoom.
	ra_c = (rra[0] + rra[1]) * 0.5;
	dec_c = (dec[0] + dec[1]) * 0.5;
	ra_w = (rra[1] - rra[0]) * 0.5;
	dec_h = (dec[1] - dec[0]) * 0.5;

	// Zoom down to 1 arc-second radius.
	while (dec_h >= (1. / 3600.)) {
		std::cout << "Frame " << fn << ", radius " << (int)floor(dec_h * 3600. + 0.5) << " arcseconds...\n";
		sprintf(fname, "frame%06d.png", fn);

		rra[0] = ra_c - ra_w;
		rra[1] = ra_c + ra_w;
		dec[0] = dec_c - dec_h;
		dec[1] = dec_c + dec_h;

		render();
		to_file(fname);

		ra_w *= 0.99;
		dec_h *= 0.99;
		++fn;
	}
}

const char* simbadasis_name(const std::string& name) {
	const char* w = name.c_str();
	if (ro.simbad_asis)
		return w;

	if (!strncmp(w, "Cl* ", 4))
		return w + 4;
	if (!strncmp(w, "Cl ", 3))
		return w + 3;
	if (!strncmp(w, "V* ", 3))
		return w + 3;
	if (!strncmp(w, "EM* ", 4))
		return w + 4;
	if (!strncmp(w, "** ", 3))
		return w + 3;
	if (!strncmp(w, "* ", 2))
		return w + 2;
	if (!strncmp(w, "NAME ", 5))
		return w + 5;

	return w;
}

int dist_label_val(double dist) {
	int ret = (int)floor(dist + 0.5);
	int step = 1;
	if (ret > 250)
		step = 5;
	if (ret > 1000)
		step = 10;
	if (ret > 4000)
		step = 20;
	if (ret > 10000)
		step = 50;
	if (ret > 20000)
		step = 100;
	if (ret > 100000)
		step = 1000;
	ret = (ret / step) * step;
	return ret;
}

char __specbuf[256];
const char* spectral_type_label(const char* st) {
	char* w;
	strcpy(__specbuf, st);
	w = __specbuf + strlen(__specbuf) - 1;
	while (w > __specbuf && !isspace(*w))
		--w;
	if (w > __specbuf)
		*w = '\000';
	return __specbuf;
}

void Map::render_name(const std::string& name, const std::string& spectral_type, int x, int y, u32 rad, double dist) {
	u32 type_size = font_size_for_height_cached(swansea, rad);
	if (type_size < 8)
		type_size = 8;
	SBitmap* render = swansea.render_cstr(simbadasis_name(name), type_size, false, nullptr);
	SBitmap* renderd = nullptr;
	if (ro.show_distance && dist > 0.) {
		char dl[512];
		if (spectral_type.empty())
			sprintf(dl, "%d ly", dist_label_val(dist));
		else
			sprintf(dl, "%d ly %s", dist_label_val(dist), spectral_type_label(spectral_type.c_str()));
		renderd = swansea.render_cstr(dl, type_size, false, nullptr);
	}
	int x1, x2, y1, y2;
	x1 = x + rad + 3;
	x2 = x1 + render->width() - 1;
	y1 = y - render->height() / 2;
	y2 = y1 + render->height() - 1;
	if (ro.show_distance && renderd != nullptr) {
		x2 = std::max(x2, x1 + (int)renderd->width() - 1);
		y2 += renderd->height() + 1;
	}
	if (isect->test_rect(x1, x2, y1, y2)) {
		Font::blit(render, map, x1, y1, ro.desig_color);
		if (renderd != nullptr)
			Font::blit(renderd, map, x1, y1 + render->height() + 2, ro.desig_color);
		isect->add_rect(x1, x2, y1, y2);
	}
	delete render;
	if (renderd != nullptr)
		delete renderd;
}

double Map::step_value(double val1, double val2) {
	double ret;
	SORT2(val1, val2, double);
	val1 = val2 - val1;

	ret = 1.0 / 3600.;
	if (val1 > (25.0 / 3600.))
		ret = 5.0 / 3600.;
	if (val1 > (75.0 / 3600.))
		ret = 10. / 3600.;
	if (val1 > (150.0 / 3600.))
		ret = 30. / 3600.;
	if (val1 > 5. / 60.)
		ret = 1.0 / 60.;
	if (val1 > 15. / 60.)
		ret = 2.0 / 60.;
	if (val1 > 30. / 60.)
		ret = 5.0 / 60.;
	if (val1 > 1.)
		ret = 10.0 / 60.;
	if (val1 > 1.75)
		ret = 15.0 / 60.;
	if (val1 > 2.5)
		ret = 30.0 / 60.;
	if (val1 > 6.0)
		ret = 1.;
	if (val1 > 12.0)
		ret = 2.;
	if (val1 > 18.0)
		ret = 3.;
	if (val1 > 24.0)
		ret = 4;
	if (val1 > 30.0)
		ret = 5.;
	if (val1 > 60.0)
		ret = 10.;

	return ret;
}

void Map::to_file(const char* fname) {
	if (!ro.no_map) {
		if (ro.downsample_width > 0) {
			SBitmap* mapds = map->resize(ro.downsample_width, 0);
			mapds->save_bmp(fname);
			delete mapds;
		} else {
			map->save_bmp(fname);
		}
	}
	if (ro.out_cone) {
		out_cone();
		out_cluster_csv();
	}
}


void Map::out_cluster_csv(void) {
	std::ofstream o;
	o.open("clusters.csv");
	o << "SOURCE_ID,X,Y,X_MILLION,Y_MILLION,DX_1000,DY_1000,COMMON_NAME\n";
	for (const auto& md : map_data) {
		int x, y, x2, y2;
		if (md.pmra == 0. && md.pmdec == 0.)
			continue;
		ra_dec_to_xy(md.ra, md.dec, x, y);
		o << md.s_id << "," << x << "," << y << ",";
		ra_dec_to_xy(md.ra + md.pmra / 3.6, md.dec + md.pmdec / 3.6, x2, y2);
		o << x2 << "," << y2 << ",";
		o << double(x2 - x) / 1000. << "," << double(y2 - y) / 1000. << ",";
		if (!ro.use_tic && !strncmp(md.desig[0].c_str(), "TIC ", 4))
			o << "\n";
		else
			o << md.desig[0] << "\n";
	}
	o.close();
	o.clear();
}

void Map::out_cone(void) {
	std::ofstream o;
	u32 cs = 0;
	o.open("cone_map.csv");
	o << "SOURCE_ID,RA,DEC,X,Y,Z,APP_MAG,DISTANCE,PX/PXE,SPECTRAL_TYPE,CATEGORY,COMMON_NAME\n";
	o.precision(std::numeric_limits<double>::max_digits10);
	for (const auto& md : map_data) {
		int x, y;
		ra_dec_to_xy(md.ra, md.dec, x, y);
		if (x < 0 || x >= w)
			continue;
		if (y < 0 || y >= h)
			continue;
		++cs;
		o << md.s_id << "," << md.ra << "," << md.dec << "," << x << "," << y << ",0," << md.MagnitudeWeighted() << "," << md.dist << ",20," << md.category << "," << md.spectral_type << ",";
		if (!ro.use_tic && !strncmp(md.desig[0].c_str(), "TIC ", 4))
			o << "\n";
		else
			o << md.desig[0] << "\n";
	}
	o.close();
	o.clear();
	std::cout << "*** " << cs << " sources within plot area output to cone_map.csv.\n";
}

void Map::ra_dec_to_xy(double ra_in, double dec_in, int& x_out, int& y_out) const {
	ra_in = (ra_in - rra[0]) / (rra[1] - rra[0]);
	dec_in = (dec_in - dec[0]) / (dec[1] - dec[0]);
	ra_in *= double(w);
	dec_in *= double(h);
	ra_in = double(w) - ra_in;
	dec_in = double(h) - dec_in;
	ra_in = floor(ra_in + 0.5);
	dec_in = floor(dec_in + 0.5);
	x_out = (int)ra_in;
	y_out = (int)dec_in;
}

void read_two_doubles(double& d1, double& d2) {
	std::string inp;
	const char* w;
	forever {
		do {
			getline(std::cin, inp);
		} while (inp.length() == 0);
		w = inp.c_str();
		d1 = atof(w);
		w = strchr(w, ',');
		if (nullptr == w) {
			std::cout << "Please enter two values separated by a comma: ";
			continue;
		}
		++w;
		d2 = atof(w);
		break;
	}
}

void read_one(char* w, double& val) {
	char* w2;
	double h = 0., m = 0., s = 0.;
	bool sign = false;

	if (strchr(w, 'h') || strchr(w, 'H')) {
		// Hours/minutes/seconds format
		__strlwr(w);
		w2 = strchr(w, 'h');
		*w2 = '\000';
		h = atof(w);
		w = w2 + 1;
		w2 = strchr(w, 'm');
		if (w2) {
			*w2 = '\000';
			m = atof(w);
			w = w2 + 1;
		}
		w2 = strchr(w, 's');
		if (w2) {
			s = atof(w);
		}
LHoursConvert:	if (h < 0.) {
			sign = true;
			h = -h;
		}
		val = (h * 15.0);
		val += (m / 4.0);
		val += (s / 240.0);	
		if (sign)
			val = -val;
		return;
	}
	if (strchr(w, ':')) {
		// Alternate HMS format.
		w2 = strchr(w, ':');
		*w2 = '\000';
		h = atof(w);
		w = w2 + 1;
		m = atof(w);
		w2 = strchr(w, ':');
		if (w2)
			s = atof(w2 + 1);
		goto LHoursConvert;
	}
	if (strchr(w, 'd') || strchr(w, 'D')) {
		// Degrees/minutes/seconds format.
		__strlwr(w);
		w2 = strchr(w, 'd');
		*w2 = '\000';
		h = atof(w);
		w = w2 + 1;
		w2 = strchr(w, 'm');
		if (w2) {
			*w2 = '\000';
			m = atof(w);
			w = w2 + 1;
		}
		w2 = strchr(w, 's');
		if (w2)
			s = atof(w);
		if (h < 0) {
			h = -h;
			sign = true;
		}
		val = h + (m / 60.) + (s / 3600.);
		if (sign)
			val = -val;
		return;
	}

	// If we're here, assume it's just decimal degrees.
	val = atof(w);
	return;
}

void read_range(char* w, double& v1, double& v2) {
	char* w2;
	w2 = strchr(w, ',');
	if (is_null(w2)) {
		read_one(w, v1);
		v2 = v1;
		return;
	}
	*w2 = '\000';
	read_one(w, v1);
	read_one(w2 + 1, v2);
}

bool object_coord(const char* name, double& ra, double& dec) {
	char line[1024];
	char lwr[256];
	std::ifstream i;
	bool ret = false;

	lwr[255] = '\000';
	strncpy(lwr, name, 255);
	__strlwr(lwr);
	i.open("xmatch/locations.csv");
	forever {
		i.getline(line, 1024);
		if (i.eof())
			break;
		__strlwr(line);
		if (!strncmp(line, lwr, strlen(lwr))) {
			char* w;
			w = strchr(line, ',') + 1;
			ra = atof(w);
			w = strchr(w, ',') + 1;
			dec = atof(w);
			ret = true;
			break;
		}
	}
	i.close();
	i.clear();
	return ret;
}

#define	valid_coord(x)	((x) >= -90. && (x) <= 360.)

int main(int argc, char* argv[]) {
	double ra1 = -999., ra2 = -999., radius = 0., dec1 = -999., dec2 = -999.;
	u32 w = 0, h = 0;
	bool dozoom = false;

	std::cout << "GaiaMap -- a star chart plotter by Chris Street. Copyright (c) 2022 Chris Street.\n";
	for (int e = 1; e < argc; ++e) {
		if (!strcmp(argv[e], "-help") || !strcmp(argv[e], "--help") || !strcmp(argv[e], "-?") || !strcmp(argv[e], "/?")) {
			std::cout << "gaiamap command line options:\n";
			std::cout << "\t-close [val]\tSet the maximum source distance (sources at least this close.)\n";
			std::cout << "\t-dec [v(s)] \tSpecify declination range, or a central DEC with -radius.\n";
			std::cout << "\t-decadj [v] \tUse this value in cos(DEC) adjustment (for tiling non-rectangular proj.)\n";
			std::cout << "\t-dim        \tDim the color of more distant sources.\n";
			std::cout << "\t-dist       \tDraw distance estimate (in light years) for sources with data.\n";
			std::cout << "\t-far [val]  \tSet the minimum source distance (sources at least this far.)\n";
			std::cout << "\t-fixra      \tInstead of adjusting RA for aspect ratio, adjust width/height instead.\n";
			std::cout << "\t-h          \tSet the height, in pixels, of the rendered map.\n";
			std::cout << "\t-help       \tThis message.\n";
			std::cout << "\t-mag [val]  \tNo sources dimmer than this magnitude will be drawn.\n";
			std::cout << "\t-neg        \tDraw the chart in negative/printable colors (black on white).\n";
			std::cout << "\t-nocone     \tSuppress the output of the cone_map.csv file.\n";
			std::cout << "\t-nomap      \tNo map is output; the cone file will be output.\n";
			std::cout << "\t-nonames    \tNo source designations will be drawn.\n";
			std::cout << "\t-object [n] \tCenter the chart on the named object.\n";
			std::cout << "\t-pm         \tDraw proper motion arrows for sources with PM data available.\n";
			std::cout << "\t-pmy [val]  \tDraw this number of years of proper motion on the chart (default 1000)\n";
			std::cout << "\t-ra [val(s)]\tSpecify right ascension range, or a central RA with -radius.\n";
			std::cout << "\t-radius [v] \tSpecify the area around point RA/DEC to plot, in arc-minutes.\n";
			std::cout << "\t-rect       \tUse a purely rectagonal projection.\n";
			std::cout << "\t-short      \tRender only the shortest designation for each source.\n";
			std::cout << "\t-simbad     \tRender SIMBAD designations as-is.\n";
			std::cout << "\t-small      \tDraw smaller stars.\n";
			std::cout << "\t-tic        \tTurn off TESS Input Catalog designations.\n";
			std::cout << "\t-w          \tSet the width, in pixels, of the rendered map. (Also -width.)\n";
			std::cout << "\t-zoom       \tCreate a succession of rendered frames zooming into this location.\n";
			std::cout << "Right ascension/declination can be in decimal degree format, or in hms/dms format.\n";
			return(0);
		}
		if (!strcmp(argv[e], "-mag")) {
			// Limiting magnitude
			if (e + 1 == argc)
				break;
			++e;
			ro.mag_limit = atof(argv[e]);
			std::cout << ">>> Option: Limiting magnitude set to [" << ro.mag_limit << "].\n";
			continue;
		}
		if (!strcmp(argv[e], "-nonames")) {
			// Turn off star designations
			ro.show_desig = false;
			std::cout << ">>> Option: rendering source designations is turned off.\n";
			continue;
		}
		if (!strcmp(argv[e], "-close")) {
			// Only sources at least this close
			if (e + 1 == argc)
				break;
			++e;
			ro.close_limit = atof(argv[e]);
			std::cout << ">>> Option: Maximum source distance set to [" << ro.close_limit << " ly].\n";
			continue;
		}
		if (!strcmp(argv[e], "-far")) {
			// Only sources at least this close
			if (e + 1 == argc)
				break;
			++e;
			ro.far_limit = atof(argv[e]);
			std::cout << ">>> Option: Minimum source distance set to [" << ro.far_limit << " ly].\n";
			continue;
		}
		if (!strcmp(argv[e], "-short")) {
			// Use short designations.
			ro.short_desig = true;
			std::cout << ">>> Option: Only the shortest designation will be drawn for each source.\n";
			continue;
		}
		if (!strcmp(argv[e], "-nomap")) {
			// Nix the map render, only output the cone file.
			ro.no_map = true;
			ro.out_cone = true;
			std::cout << ">>> Option: No map will be rendered; the cone file will still be output.\n";
			continue;
		}
		if (!strcmp(argv[e], "-nocone")) {
			// No cone_map.csv output.
			ro.out_cone = false;
			std::cout << ">>> Option: The cone file will not be output.\n";
			continue;
		}
		if (!strcmp(argv[e], "-ra")) {
			// Read the right ascension from the command line. Permits specifying in degrees or in h/m/s format.
			if (e + 1 == argc)
				break;
			++e;
			read_range(argv[e], ra1, ra2);
			if (ra1 == ra2)
				std::cout << ">>> Option: Read a single right ascension value [" << ra1 << "].\n";
			else
				std::cout << ">>> Option: Read a right ascension range [" << ra1 << "," << ra2 << "].\n";
			continue;
		}
		if (!strcmp(argv[e], "-dec")) {
			// Read the declination from the command line. This also supports h/m/s format, though that's seldom used.
			if (e + 1 == argc)
				break;
			++e;
			read_range(argv[e], dec1, dec2);
			if (dec1 == dec2)
				std::cout << ">>> Option: Read a single declination value [" << dec1 << "].\n";
			else
				std::cout << ">>> Option: Read a declination range [" << dec1 << "," << dec2 << "].\n";
			continue;
		}
		if (!strcmp(argv[e], "-radius") || !strcmp(argv[e], "-r")) {
			// Read the radius (in minutes of arc) to use to determine the plotting range.
			if (e + 1 == argc)
				break;
			++e;
			radius = atof(argv[e]);
			std::cout << ">>> Option: Radius of " << radius << " arcminutes set for use with lone RA/DEC value.\n";
			continue;
		}
		if (!strcmp(argv[e], "-w") || !strcmp(argv[e], "-width")) {
			// Set the chart area width.
			if (e + 1 == argc)
				break;
			++e;
			w = (u32) atof(argv[e]);
			std::cout << ">>> Option: Chart area pixel width set to " << w << ".\n";
			continue;
		}
		if (!strcmp(argv[e], "-h") || !strcmp(argv[e], "-height")) {
			// Set the chart area height.
			if (e + 1 == argc)
				break;
			++e;
			h = (u32) atof(argv[e]);
			std::cout << ">>> Option: Chart area pixel height set to " << h << ".\n";
			continue;
		}
		if (!strcmp(argv[e], "-dim")) {
			// Dim the color of distant sources.
			ro.dim_far = true;
			std::cout << ">>> Option: Distant sources will have their colors dimmed.\n";
			continue;
		}
		if (!strcmp(argv[e], "-neg")) {
			// Negative colors (for printing).
			ro.neg_printable = true;
			ro.bg_color = C_WHITE;
			ro.grid_color = C_GRAY;
			ro.desig_color = C_GRAY;
			ro.coord_color = C_GRAY;
			std::cout << ">>> Option: Negative (printable) colors are on.\n";
			continue;
		}
		if (!strcmp(argv[e], "-small")) {
			// Use the smaller stars.
			ro.small_stars = true;
			std::cout << ">>> Option: Will render smaller stars.\n";
			continue;
		}
		if (!strcmp(argv[e], "-object") || !strcmp(argv[e], "-obj") || !strcmp(argv[e], "-o")) {
			// Get the coordinates from a named object's location.
			if (e + 1 == argc)
				break;
			++e;
			if (object_coord(argv[e], ra1, dec1)) {
				std::cout << ">>> Option: Location set to RA = " << ra1 << ", DEC = " << dec1 << " for object \"" << argv[e] << "\".\n";
				ra2 = ra1;
				dec2 = dec1;
			} else {
				std::cout << "*** NOTICE: Object \"" << argv[e] << "\" not found in locations list.\n";
			}
			continue;
		}
		if (!strcmp(argv[e], "-simbad")) {
			// Render SIMBAD designations as-is.
			ro.simbad_asis = true;
			std::cout << ">>> Option: SIMBAD designations will be rendered as-is, including internal prefixes.\n";
			continue;
		}
		if (!strcmp(argv[e], "-rect")) {
			// Use a purely rectagonal projection, with no cos(DEC) adjustment.
			ro.rect_proj = true;
			std::cout << ">>> Option: Will use a rectagonal projection.\n";
			continue;
		}
		if (!strcmp(argv[e], "-fixra")) {
			// Fix the right ascension range; adjust width or height to maintain proper aspect ratio.
			ro.fix_ra = true;
			std::cout << ">>> Option: RA range will be fixed; unconstrained width or height will be adjusted.\n";
			continue;
		}
		if (!strcmp(argv[e], "-tic")) {
			// Ignore TIC designations.
			ro.use_tic = false;
			std::cout << ">>> Option: Will not render or output to cone TESS Input Catalog designations.\n";
			continue;
		}
		if (!strcmp(argv[e], "-decadj")) {
			// Set declination value to use for cos(DEC) adjustment.
			if (e + 1 == argc)
				break;
			++e;
			ro.use_dec = atof(argv[e]);
			if (ro.use_dec < 90. && ro.use_dec > -90.)
				std::cout << ">>> Option: Will use the DEC value " << ro.use_dec << " degrees in making cos(DEC) adjustment.\n";
			else
				std::cout << "*** Warning: invalid -decadj value " << ro.use_dec << " will be ignored; should be in (-90, 90).\n";
			continue;
		}
		if (!strcmp(argv[e], "-pm")) {
			// Draw proper motion arrows on the chart.
			ro.draw_pm = true;
			std::cout << ">>> Option: Will draw proper-motion arrows for sources with data on the chart.\n";
			continue;
		}
		if (!strcmp(argv[e], "-pmy")) {
			// Draw proper motion for a specified number of years on the chart.
			if (e + 1 == argc)
				break;
			++e;
			ro.draw_pm = true;
			ro.pm_years = atof(argv[e]);
			std::cout << ">>> Option: Will draw " << ro.pm_years << " years' worth of proper motion on the chart.\n";
			continue;
		}
		if (!strcmp(argv[e], "-dist")) {
			ro.show_distance = true;
			std::cout << ">>> Option: Will render distance estimates and spectral types on the chart.\n";
			continue;
		}
		if (!strcmp(argv[e], "-dark")) {
			ro.bg_color = C_BLACK;
			std::cout << ">>> Option: background color set to black.\n";
			continue;
		}
		if (!strcmp(argv[e], "-nogrid")) {
			ro.draw_grid = false;
			std::cout << ">>> Option: will not draw coordinate grid.\n";
			continue;
		}
		if (!strcmp(argv[e], "-downsample")) {
			if (e + 1 == argc)
				break;
			++e;
			ro.downsample_width = atoi(argv[e]);
			std::cout << ">>> Option: will downsample to width " << ro.downsample_width << ".\n";
			continue;
		}
		if (!strcmp(argv[e], "-render")) {
			ro.draw_grid = false;
			ro.show_desig = false;
			ro.draw_pm = false;
			ro.show_distance = false;
			ro.bg_color = C_BLACK;
			std::cout << ">>> Option: In render mode, -dark -nonames -nogrid -nopm -nodist set.\n";
			continue;
		}
		if (!strcmp(argv[e], "-zoom")) {
			dozoom = true;
			std::cout << ">>> Option: will render multiple frames as a zoom in.\n";
		}
	}

	if (valid_coord(ra1) && valid_coord(dec1) && radius > 0.) {
		std::cout << "SIMBAD query for objects in this region:\n";
		std::cout << "simbad.u-strasbg.fr/simbad/sim-coo?Coord=" << ra1 << "+" << dec1 << "+&Radius=" << radius << "&Radius.unit=arcmin&submit=submit+query&output.format=ASCII\n";
	}
	if (valid_coord(ra1) && ra1 == ra2 && radius > 0.) {
		ra1 -= (radius / 60.);
		ra2 += (radius / 60.);
		if (ra1 < 0.)
			ra1 = 0.;
		if (ra2 > 360.)
			ra2 = 360.;
		std::cout << "Right ascension range calculated: [" << ra1 << "," << ra2 << "].\n";
	}
	if (valid_coord(dec1) && dec1 == dec2 && radius > 0.) {
		dec1 -= (radius / 60.);
		dec2 += (radius / 60.);
		if (dec1 < -90.)
			dec1 = -90.;
		if (dec2 > 90.)
			dec2 = 90.;
		std::cout << "Declination range calculated: [" << dec1 << "," << dec2 << "].\n";
	}
	if (ra1 == ra2 || !valid_coord(ra1)) {
		std::cout << "Enter right ascension range, separated by a comma: ";
		read_two_doubles(ra1, ra2);
	}
	if (dec1 == dec2 || !valid_coord(dec1)) {
		std::cout << "Enter declination range, separated by a comma    : ";
		read_two_doubles(dec1, dec2);
	}

	if (0 == w && 0 == h) {
		char buf[16];
		std::cout << "[You can set 0 for width or height to let the aspect ratio determine that value.]\n";
		std::cout << "Enter the width of the map, in pixels : ";
		std::cin >> w;
		std::cout << "Enter the height of the map, in pixels: ";
		std::cin >> h;
		// eat the new line
		std::cin.getline(buf, 16);
	}

	Map map(ra1, ra2, dec1, dec2, w, h);
	if (dozoom) {
		map.render_zoom_in();
	} else {
		map.render();
		map.to_file("output.png");
	}

	return 0;
}


/*** end gaiamap.cpp ***/
