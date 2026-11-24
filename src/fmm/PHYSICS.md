# The physics of the Fast Multipole Method

A qualitative picture, with no equations. The central object is the gravitational potential: the field that all the masses together create in space. A body moves because of the slope of that field where it sits. The whole method is a way of knowing the potential at every body without asking every other body one at a time.

For the stages and the names used in the code, see `README.md` in this folder.

## Why the direct way wastes effort

Every body adds its own dimple to the potential, and the dimples add. To find the potential at one body, you add up the dimples of every other body. Do that for every body, and the work grows with the square of the number of bodies.

Most of that work is wasted. A far-off group of bodies produces nearly the same field at your location as its neighbours do, so you keep adding up almost the same thing.

## Distance washes out detail

Look at a distant group of masses. From far away you cannot tell the individual bodies apart. The group looks like one lump with a certain total mass. At closer range you begin to see that the mass is lopsided, then elongated, then more complicated.

Each finer feature of the arrangement fades with distance faster than the one before it. So a short list of the group's features describes its field well from far off:

- its total mass,
- then how lopsided it is,
- then how elongated it is,
- and so on.

The farther away you are, the shorter the list needs to be.

## The same fact from the receiving end

Now stand inside a small region and ask what potential the distant masses create there. Over a small region it varies smoothly:

- at lowest order it is a uniform pull, the same everywhere in the region,
- next is a tidal field, a gentle stretching that differs slightly from place to place,
- then finer variations.

The Earth feels the Moon this way. The Moon's pull on the Earth is almost uniform, and the tidal stretching that remains is what raises the ocean tides. Every body in the region feels the same ambient field and just reads off its own local value. Nobody needs to know where the Moon's individual rocks are.

So there are two short descriptions of a potential:

| description | says | valid |
|---|---|---|
| outgoing | what a group of masses looks like from far away | outside the group |
| ambient | what the field from far-away masses looks like inside a region | inside the region |

## What the algorithm does

For any region, the rest of the world falls into two zones:

- **Near neighbours.** Their dimples are sharp and detailed here, so no summary is accurate. These are treated exactly, body by body, with ordinary Newtonian gravity.
- **The far world.** Its field is smooth here, so it can be described by the short ambient description.

"Near" and "far" are measured against the size of the region, which makes the picture multi-scale. A large region sees only the very distant world as smooth. A smaller region inside it sees *more* of the world as smooth, because it is small compared with the distance. But everything that was already smooth over the large region is also smooth over the small one, so the small region *inherits* that ambient field from the larger one and only has to add what has newly become smooth at its own scale.

The potential at a body is therefore the sum of:

- the ambient fields inherited and added at every scale, from the largest regions down to its own,
- the exact contributions of its immediate neighbours.

## How the information flows

1. **Up.** Describe each small group of bodies by how it looks from afar. Bigger groups are described by combining their parts, because fields superpose. The reference point has to move, though: a lump that sits away from a new reference point looks lopsided to someone standing there, even if it looked perfectly round from its own center.
2. **Across.** Each region asks of each distant group that is now smooth over it: what ambient field do you create here? This converts the outgoing description into the ambient one.
3. **Down.** Each large region hands its ambient field to the regions inside it. It is the same field, now measured from a different reference point.
4. **Finally.** Each body feels the ambient field at its own position, plus the exact pull of its near neighbours.

## Why it is fast

Nobody talks to everyone. A region only ever exchanges descriptions with a bounded number of partners at its own scale, and everything farther arrives through its larger surroundings. Information about the field moves up, across and down, not body to body.

## What is approximate

Nothing is averaged out, and no body is ignored. The split into near and far is exact, since the potential is a straight sum of all the dimples. The only approximation is truncation: both descriptions are cut off at some finest detail, and features finer than that are discarded.

Those features fade quickly with distance, which is why the error shrinks as you keep more terms, and why it stays small as long as a region is small compared with how far away the sources are.

## Where each part lives in the code

| flow | in the code |
|---|---|
| up | `accumulate_mass_distribution` (P2M, then M2M) |
| across | `assemble_potential` (M2L for the far world, P2P for the near neighbours) |
| down | the downward pass (L2L, then L2P) |
