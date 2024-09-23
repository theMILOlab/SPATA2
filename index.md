# SPATA2 Version 3

We are excited to announce the release of SPATA2 v3.0.0, coinciding with the official publication:

Kueckelhaus, J., Frerich, S., Kada-Benotmane, J. et al. Inferring histology-associated gene expression gradients in spatial transcriptomic studies. Nat Commun 15, 7280 (2024). https://doi.org/10.1038/s41467-024-50904-x

You can install the new version by following the brief tutorial available under the Installation tab.

While the homepage is still under construction, the rest of the website, particularly the tutorials, has already been updated. Notably, we've refined the spatial gradient screening approach during the revision process. In addition, SPATA2 v3.0.0 offers numerous new features and expanded support for platforms beyond Visium.

As outlined in the newslog, SPATA2 introduces significant changes to the architecture of SPATA2 objects. The `updateSpataObject()` function should help ease this transition. Use it like this: 

`object_new <- updateSpataObject(object_old)`

However, if you encounter any issues, please don’t hesitate to open an issue, and we will do our best to assist you in transferring your analysis progress from v2.0.4 to v3.0.0. Lastly, we apologize for any delayed responses to issues over the past few months. Our team is small, and the transition from v2.0.4 to v3.0.0 required considerable time and effort. With the official publication now complete, we are committed to responding to your questions more quickly and reliably. Thank you for your patience.

## Licences Information
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
